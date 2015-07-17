# Transforms a standoff file into SBML
# For more information see (...)
# Based on information in
# [1] Buechel, F., Wrzodek, C., Mittag, F., Draeger, A., Eichner, J., Rodriguez, N., Le Novere, N., and Zell, A. (2012). Qualitative translation of relations from biopax to sbml qual. Bioinformatics, 28(20):2648?2653.
# [2] Pyysalo, S., Ohta, T. and Ananiadou, S. (2012) Pathway Curation (PC) task at BioNLP Shared Task 2013 Task Proposal

# to run all tests
# ls tests/*.ann | xargs -I {} python2.7 st2sbml.py --file {}  --output {}.output-sbml.xml --complete-reactions

import argparse
import cgi
import copy
import logging
import os
import re
import sets
import sys

import libsbml
import parse_standoff
import uniprot_tools

# mapping of trigger (entity) to sbo term (based on 1)
STANDOFF_ENTITY_TO_SBO_MAPPING = {
    "gene" : "SBO:0000354", # informational molecule segment
    "complex" : "SBO:0000253", # non-covalent complex
    "protein" : "SBO:0000252", # polypeptide chain
    "dna" : "SBO:0000251", # deoxyribonucleic acid
    "dnaregion" : "SBO:0000251", # deoxyribonucleic acid
    "rna" : "SBO:0000250", # ribonucleic acid
    "rnaregion" : "SBO:0000250", # ribonucleic acid
    "smallmolecule" : "SBO:0000247", # simple chemical
    "simple_molecule" : "SBO:0000247", # simple chemical
    "ion" : "SBO:0000247", # simple chemical
    "drug" : "SBO:0000247" # simple chemical
}

# mapping of reaction types to sbo term 
# for reactions to whom we do not know the SBO term
# use the generic reaction SBO:0000176
BIOCHEMICAL_REACTION ="SBO:0000176"

STANDOFF_EVENT_TO_SBO_MAPPING = {
    "acetylation"  :   "SBO:0000215",
    "activation"   :   BIOCHEMICAL_REACTION,
    "association"  :   "SBO:0000297",
    "binding"  :   "SBO:0000297",
    "catabolism"   :   "GO:0009056",
    "catalysis"    :   "SBO:0000172",
    "conversion"   :   BIOCHEMICAL_REACTION,
    "deacetylation"    :   "GO:0006476",
    "degradation"  :   "SBO:0000179",
    "demethylation"    :   "GO:0006482",
    "dephosphorylation"    :   "SBO:0000330",
    "deubiquitination" :   "GO:0016579",
    "dissociation" :   "SBO:0000180",
    "gene_expression" : "GO:0010467",
    "inactivation" :   BIOCHEMICAL_REACTION,
    "localization" :   "GO:0051179",
    "methylation"  :   "SBO:0000214",
    "negative_regulation"  :   "SBO:0000169",
    "pathway"  :   "SBO:0000375", 
    "phosphorylation"  :   "SBO:0000216",
    "positive_regulation"  :   "SBO:0000170",
    "protein_catabolism"   :   "SBO:0000179", # degradation
    "regulation"   :   "SBO:0000168",
    "transcription"    :   "SBO:0000183", 
    "translation"  :   "SBO:0000184",
    "transport"    :   "SBO:0000185",
    "ubiquitination"   :   "SBO:0000224"
}

# mapping of general reaction components to sbo term 
GENERIC_REACTION_SBO_MAPPING = {
    "reactant"  :   "SBO:0000010",
    "product"  :   "SBO:0000011",
    "modifier"  :   "SBO:0000019",
    "activator"  :   "SBO:0000021",
    "inhibitor"  :   "SBO:0000020"
}

def parse_arguments( args = None):
    
    parser = argparse.ArgumentParser( description='Converts BioNLP shared task format to SBML (requires libsbml).')
    
    parser.add_argument( '--file',
                    action = "store",
                    dest = "file",
                    default = "test.ann",
                    help = "the file to convert (default is test.ann)")

    parser.add_argument( '--path',
                        action = "store",
                        dest = "path",
                        default = None,
                        help = "the path to a1, a2 or ann file (includes the prefix)")
    
    parser.add_argument( '--a1',
                        action = "store",
                        dest = "a1",
                        default = ".a1",
                        help = "the file ending for a1 files")
    
    parser.add_argument( '--a2',
                        action="store",
                        dest = "a2",
                        default = ".a2",
                        help = "the file ending for a2 files")
    
    parser.add_argument( '--ann',
                        action="store",
                        dest = "ann",
                        default = ".ann",
                        help = "the file ending for a2 files")
    
    parser.add_argument( '--use-ann',
                        action="store_true",
                        dest = "use_ann",
                        default = True,
                        help="whether to use ann or a1/a2")
    
    parser.add_argument( '--output',
                        action = "store",
                        dest = "output",
                        default = "output-sbml.xml",
                        help="the path to output file")

    parser.add_argument( '--remove-unconnected-species',
                        action = "store_true",
                        dest = "remove_unconnected_species",
                        default = False,
                        help = "whether to keep unconnected species (i.e. species not used in reactions)")

    parser.add_argument( '--complete-reactions',
                        action = "store_true",
                        dest = "complete_reactions",
                        default = False,
                        help = "whether to complete reactions by adding reactants and products.")

    parser.add_argument( '--remove-unconnected-reactions',
                        action = "store_true",
                        dest = "remove_unconnected_reactions",
                        default = False,
                        help = "whether to keep reactions without reactants, products or modifiers.")

    parser.add_argument( '--remove-empty-compartments',
                        action = "store_true",
                        dest = "remove_empty_compartments",
                        default = False,
                        help="whether to keep compartments without species in them (default is keep them).")
    
    parser.add_argument( '--default-compartment-name',
                    action ="store",
                    dest = "default_compartment_name",
                    default = "default",
                    help = "what to use as the default compartment (default is \"default\")")
    
    parser.add_argument( '--use-uniprot',
                    action = "store_true",
                    dest = "use_uniprot",
                    default = False,
                    help = "whether to add uniprot information (default is no). requires an internet connection.")

    parser.add_argument( '--log-level',
                    action = "store",
                    dest = "log_level",
                    default = "INFO",
                    help = "the log level for messages (default is INFO, can be DEBUG, INFO, WARNING, ERROR, NONE)")

    parser.add_argument( '--uniprot-pickle', action="store",
                    dest = "uniprot_pickle",
                    default = None,
                    help = "whether to get uniprot information from a pickle file and where (default is None, But you can use\"uniprot.pickle\" which is generated by uniprot_tools.py)")
    return parser.parse_args( args = args);

DEFAULT_ARGUMENTS = parse_arguments( [])

def get_path( arguments):
    if arguments.path != None:
        return arguments.path
    else:
        return arguments.file

def check( value, message, arguments = DEFAULT_ARGUMENTS):
    """If 'value' is None, prints an error message constructed using
    'message' and then exits with status code 1.  If 'value' is an integer,
    it assumes it is a libSBML return status code.  If the code value is
    LIBSBML_OPERATION_SUCCESS, returns without further action; if it is not,
    prints an error message constructed using 'message' along with text from
    libSBML explaining the meaning of the code, and exits with status code 1. """
    if value == None:
        err_msg = "Error processing {0} libSBML returned a null value trying to ".format( get_path( arguments), message)
        raise SystemExit( err_msg)
    elif type( value) is int:
        if value == libsbml.LIBSBML_OPERATION_SUCCESS:
            return
        else:
            err_msg = 'Error processing {0} tying to {1} libSBML returned error code {2} : "{3}"'.format( get_path( arguments), message, value, libsbml.OperationReturnValue_toString(value).strip())
            raise SystemExit(err_msg)
    else:
        return

def add_cvterm( term, arguments = DEFAULT_ARGUMENTS):
    "Adds a controlled vocabulary term"
    controlledVocab = libsbml.CVTerm(libsbml.BIOLOGICAL_QUALIFIER);
    check( controlledVocab.setBiologicalQualifierType(libsbml.BQB_IS), "addSBO term-pt1", arguments = arguments);
    check( controlledVocab.addResource(term),"addSBO term-pt2", arguments = arguments);
    return controlledVocab

def add_note( note, species, arguments = DEFAULT_ARGUMENTS):
    """ Adds a note to species (wraps the note in <p></p> and escapes the text) """
    check( species.appendNotes( "<p xmlns=\"http://www.w3.org/1999/xhtml\">{0}</p>".format( cgi.escape( note))), 'append notes', arguments = arguments);
    return species;

def add_species( trigger, model, id = None, name = None, compartment = "default", arguments = DEFAULT_ARGUMENTS):
    "Creates a species and adds it to model"
    
    species = None;
    
    species = model.createSpecies()
    check( species, 'create species', arguments = arguments);
    
    if id == None:
        check( species.setId( trigger.id), 'set species id', arguments = arguments);
        check( species.setMetaId( "metaid_0000" + trigger.id), 'set meta ID', arguments = arguments);
    else:
        check( species.setId( id), 'set species id', arguments = arguments);
        check( species.setMetaId( "metaid_0000" + id), 'set meta ID', arguments = arguments);
    
    if name == None:
        check( species.setName( trigger.text), 'set species name', arguments = arguments);
    else:
        check( species.setName( name), 'set species name', arguments = arguments);
    
    check( species.setCompartment( compartment), 'set species compartment', arguments = arguments);
    
    if trigger and STANDOFF_ENTITY_TO_SBO_MAPPING.get( trigger.type_lower):
        sbo_term = STANDOFF_ENTITY_TO_SBO_MAPPING.get( trigger.type_lower);
        check( species.setSBOTerm( sbo_term), 'set sbo terms', arguments = arguments);
        
    if not trigger is None: 
        add_note( "trigger: {0}".format( trigger.text), species, arguments = arguments);
        add_note( "start/end: {0}/{1}".format( trigger.start, trigger.end), species, arguments = arguments);
        
    return species

def modify_species( species, model, compartment = "default", arguments = DEFAULT_ARGUMENTS):
    species_ref = model.getSpecies( species.id)
    check( species_ref, 'modify species compartment', arguments = arguments);
    species_ref.setCompartment( compartment);
    check( species_ref.setCompartment( compartment), 'modify species compartment', arguments = arguments);

def add_reaction( event, model, arguments = DEFAULT_ARGUMENTS):
    reaction = model.createReaction();
    check( reaction, 'create reaction', arguments = arguments);
    check( reaction.setId( event.id), 'set event id', arguments = arguments);
    check( reaction.setName( event.type), 'set event name', arguments = arguments);
    check( reaction.setMetaId( "metaid_0000"+event.id), 'set meta ID', arguments = arguments);
    check( reaction.setReversible( False), 'make irreversible', arguments = arguments);
    if event and STANDOFF_EVENT_TO_SBO_MAPPING.get( event.type_lower) and STANDOFF_EVENT_TO_SBO_MAPPING[event.type_lower][0:3] == "SBO":
        check( reaction.setSBOTerm( STANDOFF_EVENT_TO_SBO_MAPPING[event.type_lower]), 'set sbo terms', arguments = arguments);
    check( reaction.addCVTerm( add_cvterm( STANDOFF_EVENT_TO_SBO_MAPPING[event.type_lower])), 'set controlled vocab', arguments = arguments);
    
    add_note( "trigger: {0}".format( event.trigger.text), reaction, arguments = arguments);
    add_note( "start/end: {0}/{1}".format( event.trigger.start, event.trigger.end), reaction, arguments = arguments);

    return reaction

def add_reactant( reactant_id, reaction, model, compartment = "default", arguments = DEFAULT_ARGUMENTS):
    "Adds a reactant to a reaction or makes one up (if product == None)"
    if not reactant_id:
        reaction_name = reaction.getName();
        reaction_id = reaction.getId();
    
        reactant_id = reaction_id + "_reactant_nr" + str(len( reaction.getListOfReactants()));
        reactant_name = reaction_name[0:3].lower() + "reactant";
        add_species( None, model, id = reactant_id, name = reactant_name, compartment = compartment, arguments = arguments);

    reactant_ref = reaction.createReactant()
    check( reactant_ref, 'create reactant reference', arguments = arguments);
    check( reactant_ref.setSpecies( reactant_id), 'assign reactant species', arguments = arguments);
    check( reactant_ref.setMetaId( "metaid_0000" + reactant_id), 'set meta ID', arguments = arguments);
    check( reactant_ref.addCVTerm(add_cvterm( GENERIC_REACTION_SBO_MAPPING["reactant"])), 'set controlled vocab SBO term for reactant', arguments = arguments);
    #check( reactant_ref.addCVTerm(add_cvterm( STANDOFF_ENTITY_TO_SBO_MAPPING[reactant.type])), 'set controlled vocab SBO term 2 for reactant')
    return reactant_ref

def add_product( product_id, reaction, model, compartment = "default", arguments = DEFAULT_ARGUMENTS):
    """ Adds a product or makes one up (if product == None and length of products > 0) """
    if product_id is None and len( reaction.getListOfProducts()) > 0:
        product_ref = reaction.getListOfProducts()[0]
        product_species = product_ref.getSpecies()
        return product_species
    else:
        if product_id is None:
            reaction_name = reaction.getName();
            reaction_id = reaction.getId();
            product_id = reaction_id + "_product_nr" + str(len( reaction.getListOfProducts()));
            
            if reaction.getName().lower() in ["transport", "localization"]:
                product_prefix = ""
            else:
                product_prefix = reaction_name[0:3].lower();
            
            if len( reaction.getListOfReactants()) > 0:
                reactant = model.getSpecies( reaction.getListOfReactants()[0].getSpecies())
                product_name = product_prefix + reactant.getName();
            else:
                product_name = product_prefix + "Product";
            add_species( None, model, id = product_id, name = product_name, compartment = compartment, arguments = arguments);
        
        product_ref = reaction.createProduct()
        check( product_ref, 'create product reference', arguments = arguments);
        check( product_ref.setSpecies( product_id), 'assign product species', arguments = arguments);
        check( product_ref.setMetaId( "metaid_0000" + product_id), 'set meta ID', arguments = arguments);
        check( product_ref.addCVTerm( add_cvterm( GENERIC_REACTION_SBO_MAPPING["product"])), 'set controlled vocab SBO term for product', arguments = arguments);
        # check( product_ref.addCVTerm(add_cvterm(STANDOFF_ENTITY_TO_SBO_MAPPING[product.type])), 'set controlled vocab SBO term 2 for product')
        return product_id

def add_modifier( modifier_id, reaction, model, compartment = "default", arguments = DEFAULT_ARGUMENTS):
    """ Adds a modifier or makes one up (if product == None and length of products != 1) """
    if modifier_id is None and len( reaction.getListOfModifiers()) == 1:
        return reaction.getListOfModifiers()[0].getId()
    else:
        if not modifier_id:
            reaction_name = reaction.getName();
            reaction_id = reaction.getId();
            modifier_id = reaction_id + "_modifier_nr" + str(len( reaction.getListOfModifiers()));
            modifier_name = reaction_name[0:3].lower() + "modifier";
            add_species( None, model, id = modifier_id, name = modifier_name, compartment = compartment, arguments = arguments);
        
        modifier_ref = reaction.createModifier()
        check( modifier_ref, 'create modifier reference', arguments = arguments);
        check( modifier_ref.setSpecies( modifier_id), 'assign modifier species', arguments = arguments);
        check( modifier_ref.setMetaId( "metaid_0000" + modifier_id), 'set meta ID', arguments = arguments);
        check( modifier_ref.addCVTerm(add_cvterm(GENERIC_REACTION_SBO_MAPPING["modifier"])), 'set controlled vocab SBO term for cause', arguments = arguments);
        return modifier_ref

def add_compartment( compartment_id, compartment_name, model, arguments = DEFAULT_ARGUMENTS):
    compartment_ref = model.getCompartment( compartment_id);
    if compartment_ref == None:
        compartment_ref = model.createCompartment()
        check( compartment_ref, 'create compartment', arguments = arguments);
        check( compartment_ref.setId( compartment_id), 'set compartment id', arguments = arguments);
        check( compartment_ref.setName( compartment_name), 'set compartment name', arguments = arguments);
        check( compartment_ref.setConstant(True), 'set compartment "constant"', arguments = arguments);
        check( compartment_ref.setSize(1), 'set compartment "size"', arguments = arguments);
        check( compartment_ref.setUnits('volume'), 'set compartment size units', arguments = arguments);
    return compartment_ref;

##################################################################
##### handle entities
##################################################################

def handle_entities( entities, model, compartment = "default", arguments = DEFAULT_ARGUMENTS):
    """ goes through a set of entities and adds them to the model as species,
        if they can be mapped to an SBO term """
    for trigger in entities:
        try:
            if STANDOFF_ENTITY_TO_SBO_MAPPING.get( trigger.type_lower):
                add_species( trigger, model, compartment = compartment, arguments = arguments);
        except:
            logging.getLogger( "st2sbml").error( "{0} entity {1} could not be added to model".format( get_path( arguments), trigger))
            
##################################################################
##### handle specific events
##################################################################

def handle_localization( event, model, arguments = DEFAULT_ARGUMENTS):
    """ Handles localization and transport events """
    # transportation and localization allow for AtLoc, FromLoc, ToLoc
    # and are handled specially (themes must be at least one and species)
    to_loc = None;
    from_loc = None;
    at_loc = None;

    if len(event.get_roles("toloc")) == 1:
        to_loc = event.get_roles("toloc")[0]
    elif len(event.get_roles("toloc")) != 0:
        logging.getLogger( "st2sbml").warning("{0} event {1} ToLoc not handled because more than 1 ToLoc".format( get_path( arguments), event)) 
        if not isinstance( to_loc, parse_standoff.EntityTrigger):
            logging.getLogger( "st2sbml").warning("{0} event {1} ToLoc {2} not handled because it is not an entity".format( get_path( arguments), event, to_loc))
            to_loc = None;
    if len(event.get_roles("fromloc")) == 1:
        from_loc = event.get_roles("fromloc")[0]
        if not isinstance( from_loc, parse_standoff.EntityTrigger):
            logging.getLogger( "st2sbml").warning("{0} event {1} FromLoc {2} not handled because it is not an entity".format( get_path( arguments), event, from_loc))
            from_loc = None;
    elif len(event.get_roles("fromloc")) != 0:
        logging.getLogger( "st2sbml").warning("{0} event {1} FromLoc not handled because more than 1 FromLoc".format( get_path( arguments), event))
    if len(event.get_roles("atloc")) == 1:
        at_loc = event.get_roles("atloc")[0]
        if not isinstance( at_loc, parse_standoff.EntityTrigger):
            logging.getLogger( "st2sbml").warning("{0} event {1} AtLoc {2} not handled because it is not an entity".format( get_path( arguments), event, at_loc))
            at_loc = None;
    elif len(event.get_roles("atloc")) != 0:
        logging.getLogger( "st2sbml").warning("{0} event {1} AtLoc not handled because more than 1 AtLoc".format( get_path( arguments), event))
        
    if from_loc and at_loc:
        logging.getLogger( "st2sbml").warning("{0} event {1} not handled because FromLoc and AtLoc given at the same time.".format( get_path( arguments), event))
    elif to_loc and at_loc:
        logging.getLogger( "st2sbml").warning("{0} event {1} not handled because ToLoc and AtLoc given at the same time.".format( get_path( arguments), event))
    elif from_loc is None and to_loc is None and at_loc is None:
        logging.getLogger( "st2sbml").warning("{0} event {1} not handled because no AtLoc, FromLoc or ToLoc given".format( get_path( arguments), event))
    else:
        if at_loc: # at_loc
            for theme in event.get_roles( "theme"):
                if isinstance( theme, parse_standoff.EntityTrigger):
                    add_compartment( at_loc.id, at_loc.text, model = model, arguments = arguments)
                    modify_species( theme, model, compartment = at_loc.id, arguments = arguments);
                else:
                    logging.getLogger( "st2sbml").warning("{0} event {1} theme {2} not handled because not a species".format( get_path( arguments), event, theme))
        else: # from_loc and/or to_loc
            reaction = add_reaction( event, model, arguments = arguments);
            for theme in event.get_roles("theme"):
                if not isinstance( theme, parse_standoff.EntityTrigger):
                    logging.getLogger( "st2sbml").warning("{0} event {1} not handled because {2} is not a species".format( get_path( arguments), event, theme))
                elif model.getSpecies( theme.id) is None:
                    logging.getLogger( "st2sbml").warning("{0} event {1} not handled because {2} not found in the model.".format( get_path( arguments), event, theme))
                else:
                    add_reactant( theme.id, reaction, model, arguments = arguments);
                if from_loc:
                    add_compartment( from_loc.id, from_loc.text, model = model, arguments = arguments);
                    modify_species( model.getSpecies( theme.id), model = model, compartment = from_loc.id, arguments = arguments);
                if to_loc:
                    add_compartment( to_loc.id, to_loc.text, model, arguments = arguments);
                    product_id = add_product( product_id = None, reaction = reaction, compartment = to_loc.id, model = model, arguments = arguments);
                    modify_species( model.getSpecies( product_id), model = model, compartment = to_loc.id, arguments = arguments);

def handle_regulation( event, model, arguments = DEFAULT_ARGUMENTS):
    unhandled_roles = []
    
    # we add reaction if there is a theme which is not an event
    if len( event.get_roles( "Theme")) == 0 or any( [isinstance( theme, parse_standoff.EntityTrigger) for theme in event.get_roles( "Theme")]):
        reaction = add_reaction( event, model, arguments = arguments);
        
        for theme in event.get_roles( "theme"):
            if isinstance( theme, parse_standoff.EntityTrigger):
                add_reactant( theme.id, reaction, model, arguments = arguments);
            elif theme.id in isinstance( cause, parse_standoff.Event):
                unhandled_roles.append( ( "theme", theme))
            else:
                logging.getLogger( "st2sbml").warning( "{0} event {1} theme {2} unhandled.".format( get_path( arguments), event, theme))
        for cause in event.get_roles("cause"):
            if isinstance( cause, parse_standoff.EntityTrigger):
                add_modifier( cause.id, reaction, model, arguments = arguments);
            elif isinstance( cause, parse_standoff.Event):
                unhandled_roles.append( ( "cause", cause))
            else:
                logging.getLogger( "st2sbml").warning( "{0} event {1} cause {2} unhandled.".format( get_path( arguments), event, theme))
    else:
        unhandled_roles = event.roles;

    # return unhandled roles
    if len(unhandled_roles) > 0:
        return { event.id : unhandled_roles }
    else:
        return None;

def handle_gene_expression( event, model, arguments = DEFAULT_ARGUMENTS):
    """ Handle Gene Expression, Transcription and Translation """
    # Transcription RNA from nothing- Caused by Gene
    # Translation Protein from nothing - Caused by RNA
    reaction = add_reaction( event, model, arguments = arguments);
    # for translation proteins are products (everything else is modifier)
    if event.type_lower == "translation":
        for theme in event.get_roles("theme"):
            if theme.type == "Protein":
                add_product( theme.id, reaction, model, arguments = arguments);
            else:
                add_modifier( theme.id, reaction, model, arguments = arguments);
    # for gene_expression and transcription - Rna and proteins are products
    else:
        for theme in event.get_roles("theme"):
            if theme.type_lower == "rna" or theme.type_lower == "protein":
                add_product( theme.id, reaction, model, arguments = arguments);
            else:
                add_modifier( theme.id, reaction, model, arguments = arguments);
                
##################################################################
##### HANDLE UNHANDLED EVENTS (events where some role is an event)
##################################################################

def resolve_cause_ids( cause, model, arguments = DEFAULT_ARGUMENTS):
    """ takes an entity or event (which is the cause of a regulation event and retrieves the actual cause id) """
    # cause can be an entity -> cause is the entity itself
    if isinstance( cause, parse_standoff.EntityTrigger):
        return [cause.id]
    # cause can be an event which is in the model
    # -> then the cause is actually the product of that event
    elif model.getReaction( cause.id):
        reaction = model.getReaction( cause.id)
        product_id = add_product(  product_id = None, reaction = reaction, model = model, arguments = arguments);
        return [ product_id]
    # cause can be an event which is Regulation and not in model
    # -> cause is the cause of that event
    elif  isinstance( cause, parse_standoff.Event) \
        and cause.type in ["positive_regulation", "negative_regulation", "regulation", "catalysis"] \
        and len( cause.get_roles("cause")) > 0:
        # find the causes of the cause event
        results = []
        for c in cause.get_roles("cause"):
            results.extend( resolve_cause_ids( c, model, arguments = arguments))
        return results
    elif  isinstance( cause, parse_standoff.Event) \
        and cause.type in ["positive_regulation", "negative_regulation", "regulation", "catalysis"]:
        cause_entity = add_species( None,
                                   model = model,
                                   id = cause.id + "_Cause_0",
                                   name = "Cause",
                                   arguments = arguments);
        return [ cause_entity.getId()]
    # cannot handle other causes
    else:
        return [];
    
def resolve_themes( theme, model):
    """ takes an event (which is the theme of a regulation event and retrieves it's event) """
    ## finds the theme can either be an entity or an event that exists in the model
    # cause can be an entity or an event
    if  model.getReaction( theme.id):
        return [theme]
    # theme can be a Regulation event not added to the model because
    # it consists of event roles
    elif theme.type in ["regulation", "negative_regulation", "positive_regulation", "catalysis"]:
        themes = []
        for theme in theme.get_roles("theme"):
            themes.extend( resolve_themes( theme, model))
        return themes
    else:
        return [];
    
def handle_unhandled_event( event, unhandled_roles, model, arguments = DEFAULT_ARGUMENTS):
    """ Infers real event theme and causes and add them to the model if possible """
    causes = [ role[1] for role in unhandled_roles if role[0].lower().startswith("cause") ]
    themes = [ role[1] for role in unhandled_roles if role[0].lower().startswith("theme") ]
    
    # find the actual themes of the event
    unhandled_themes = []
    for theme in themes:
        ts = resolve_themes( theme, model)
        if ts == []:
            logging.getLogger( "st2sbml").warning( "{0} event {1} theme {2} not handled (unable to resolve reaction)".format( get_path( arguments), event, theme))
        unhandled_themes.extend( ts)
    all_themes = set( event.get_roles( "theme") + unhandled_themes)
        
    # find the actual causes of the event
    unhandled_cause_ids = [];
    for cause in causes:
        cs = resolve_cause_ids( cause, model, arguments = arguments)
        if cs == []:
            logging.getLogger( "st2sbml").warning( "{0} event {1} cause {2} not handled (unable to resolve species)".format( get_path( arguments), event, cause))
        unhandled_cause_ids.extend( cs)
    all_cause_ids = set([ cause.id for cause in event.get_roles( "Cause") if isinstance( cause, parse_standoff.EntityTrigger)] + unhandled_cause_ids)

    reaction = model.getReaction( event.id);

    # the event does not exist in model (this means that themes are only events)
    if reaction is None:
        if unhandled_themes != [] and all_cause_ids != []:
            for cause_id in all_cause_ids:
                for theme in unhandled_themes:
                    theme_reaction = model.getReaction( theme.id);
                    add_modifier( cause_id, theme_reaction, model, arguments = arguments);
            logging.getLogger( "st2sbml").debug( "{0} event {1} is not a reaction but causes added to theme event(s)".format( get_path( arguments), event))
        elif unhandled_themes == [] and all_cause_ids == []:
            logging.getLogger( "st2sbml").warning( "{0} event {1} not added, no themes and no causes found".format( get_path( arguments), event))
        elif all_cause_ids == []:
            logging.getLogger( "st2sbml").warning( "{0} event {1} not added, no causes found".format( get_path( arguments), event))
        else:
            logging.getLogger( "st2sbml").warning( "{0} event {1} not added, no theme found".format( get_path( arguments), event))
    else: # the event exists in the model - some themes are events
        # handle unhandled themes
        for theme in unhandled_themes:
            if model.getReaction( theme.id): # event which is in Model
                product = add_product( product_id = None, reaction = reaction, model = model, arguments = arguments);
                add_modifier( product.getId(), model.getReaction( theme.id), model, arguments = arguments);
            else:
                logging.getLogger( "st2sbml").warning( "{0} event {1} theme {2} not handled (theme is not a species or not a reaction in the model)".format( get_path( arguments), event, theme))
        # handle causes
        for cause_id in unhandled_cause_ids:
            add_modifier( cause_id, reaction, model, arguments = arguments);

##################################################################
##### HANDLE EVENTS
##################################################################

def handle_events( events, model, arguments = DEFAULT_ARGUMENTS):
    """ adds events to the model as reactions
        events = dictionary of event_id : event
        model = an initialized sbml Model """
    unhandled_events = {}
    
    # pass 1 all events
    for event_id in events:
        event = events[event_id]
        
        if STANDOFF_EVENT_TO_SBO_MAPPING.get( event.type_lower) is None:
            # event is unknown
            logging.getLogger( "st2sbml").warning( "{0} event {1} unhandled, because unknown".format( get_path( arguments), event))
        elif event.type_lower == "pathway":
            pass # do nothing for pathways (entities have already been added)
        # handle localization (special handling)
        elif event.type_lower in [ "localization", "transport"]: 
            handle_localization( event, model, arguments = arguments);
        # handle regulation events (special handling)
        elif event.type_lower in ["regulation", "positive_regulation", "negative_regulation", "activation", "inactivation", "catalysis"]:
            unhandled_event = handle_regulation( event, model, arguments = arguments);
            if not unhandled_event is None:
                unhandled_events.update( unhandled_event);
        elif event.type_lower in ["gene_expression", "transcription", "translation"]:
            handle_gene_expression( event, model, arguments = arguments);
        # not all roles are entities
        elif not all( [ isinstance( role[1] , parse_standoff.EntityTrigger) for role in event.roles]):
            logging.getLogger( "st2sbml").warning( "{0} event {1} unhandled. Some roles are events, which is not allowed for this event type".format( get_path( arguments), event))
        # everything else: Conversion, Acetylation, Deacetylation, Demethylation, Dephosphorylation, Deubiquitination, Methylation, Phosphorylation, Ubiquitination
        else: 
            # add reaction
            reaction = add_reaction( event, model, arguments = arguments);
            
            # handle products -> add as product
            for product in event.get_roles( "product"):
                add_product( product.id, reaction, model, arguments = arguments);
            # handle themes -> add as reactants
            for theme in event.get_roles( "theme"):
                add_reactant( theme.id, reaction, model, arguments = arguments);
            # handle comp -> add as reactants
            for comp in event.get_roles( "complex"):
                add_reactant( comp.id, reaction, model, arguments = arguments);
            # handle Participant -> add as product
            for comp in event.get_roles( "participant"):
                add_reactant( comp.id, reaction, model, arguments = arguments);
            # handle causes -> add as modifiers
            for cause in event.get_roles( "cause"):
                add_modifier( cause.id, reaction, model, arguments = arguments);
            for site in event.get_roles( "site"):
                add_note( "Site: {0}".format( site.text), reaction, arguments = arguments);

            # check if there are any unhandled roles
            for unhandled_role in set([role[0] for role in event.roles]).difference(["theme", "cause", "product", "site", "participant", "complex"]):
                logging.getLogger( "st2sbml").warning( "{0} event {1} role {2} not handled, because unknown.".format( get_path( arguments), event, unhandled_role))


    # pass 2 all unhandled events
    for event_id in unhandled_events:
        event = events[event_id]
        unhandled_roles = unhandled_events[event_id]
        handle_unhandled_event( event, unhandled_roles, model, arguments = arguments);
    
                    
##################################################################
##### CLEANUP OPERATIONS
##################################################################

def cleanup_add_uniprot_information( model, entities, arguments = DEFAULT_ARGUMENTS):
    for entity in entities:
        species = model.getSpecies( entity.id)
        if not species is None:
            uniprot_data = uniprot_tools.get_uniprot_data( entity.text, uniprot_pickle = arguments.uniprot_pickle) 
            if len( uniprot_data.keys()) == 0:
                logging.getLogger( "st2sbml").warning( "{0} no uniprot information found for {1}".format( get_path( arguments), entity.text))
            else:
                if not uniprot_data.get("official_name") is None:
                    add_note( "original name: {0}".format( uniprot_data.get("official_name")), species, arguments = arguments);
                if not uniprot_data.get("gene_id") is None:
                    add_note( "gene id: {0}".format( uniprot_data.get("gene_id")), species, arguments = arguments);
                if not uniprot_data.get("gene_name") is None:
                    add_note( "gene name: {0}".format( uniprot_data.get("gene_name")), species, arguments = arguments);
                for alternate_name in uniprot_data.get("alternate_names"):
                    add_note( "alternate name: {0}".format( alternate_name), species, arguments = arguments);
                for gene_name_synonym in uniprot_data.get("gene_name_synonyms"):
                    add_note( "gene name synonyms: {0}".format( gene_name_synonym), species, arguments = arguments);
                for uniprot_id in uniprot_data["uniprot_ids"]:
                    check( species.addCVTerm( add_cvterm( "urn:miriam:uniprot:{0}".format(uniprot_id))), 'set controlled vocab', arguments = arguments);

def cleanup_remove_unconnected_species( model, arguments = DEFAULT_ARGUMENTS):
    #Remove unwanted species (the ones not involved in a any reaction)
    set_of_species = sets.Set( [s.id for s in model.getListOfSpecies()])
    set_of_species_in_reactions = sets.Set()
    for reaction in model.getListOfReactions():
        reactants = sets.Set( [s.getSpecies() for s in reaction.getListOfReactants()])
        modifiers = sets.Set( [s.getSpecies() for s in reaction.getListOfModifiers()]) 
        products = sets.Set( [s.getSpecies() for s in reaction.getListOfProducts()])
        reaction_species = reactants.union( products).union( modifiers)
        set_of_species_in_reactions = set_of_species_in_reactions.union( reaction_species)
        
    for species_id in set_of_species.difference( set_of_species_in_reactions):
        check( model.removeSpecies( species_id), 'removing unconnected species', arguments = arguments);

def cleanup_remove_unconnected_reactions( model, arguments = DEFAULT_ARGUMENTS):
    # remove reactions w/o reactants, products or modifiers
    reaction_ids = [ r.getId() for r in model.getListOfReactions()];
    for reaction_id in reaction_ids:
        reaction = model.getReaction( reaction_id)
        if len(reaction.getListOfReactants()) == 0 and len(reaction.getListOfProducts()) == 0 and len(reaction.getListOfModifiers()) == 0:
            check( model.removeReaction( reaction_id), 'removing reaction w/o reactants, products and modifiers', arguments = arguments);

def cleanup_complete_reactions( model, events, arguments = DEFAULT_ARGUMENTS):
    # complete reactions with reactants and products if none exist
    for reaction in model.getListOfReactions():
        if len( reaction.getListOfReactants()) == 0:
            event = events[reaction.getId()];
            if event is None or event.type_lower not in ["gene_expression", "transcription", "translation"]:
                add_reactant( None, reaction, model, arguments = arguments);
        if len( reaction.getListOfProducts()) == 0:
            event = events[reaction.getId()];
            if event is None or event.type_lower not in ["degradation", "catabolism", "protein_catabolism"]:
                add_product( None, reaction, model, arguments = arguments);

def cleanup_remove_empty_compartments( model, arguments = DEFAULT_ARGUMENTS):
    # remove empty compartments
    compartments_used = sets.Set([ species.getCompartment() for species in model.getListOfSpecies()])
    compartment_ids = sets.Set([ c.getId() for c in model.getListOfCompartments()]);
    compartments_not_used = compartment_ids.difference( compartments_used)
    for compartment_id in compartments_not_used:
        check( model.removeCompartment( compartment_id), 'removing compartment', arguments = arguments)
        
##################################################################
##### CREATE MODEL
##################################################################
    
def create_document( trigger, entity_trigger, events, arguments = DEFAULT_ARGUMENTS):
    """ Creating an SBML model from entities and events """
    ##### CREATE SBMLDocument
    try:
      document = libsbml.SBMLDocument( 2, 4)
    except ValueError:
        logging.getLogger( "st2sbml").error( 'Could not create SBMLDocument object')
        sys.exit(1)
    
    ##### CREATE the model 
    model = document.createModel()
    check( model, 'create model', arguments = arguments);
    
    # Create default compartment
    add_compartment( 'default', 'default', model, arguments = arguments);
    
    # add entities to the model as species
    handle_entities( entity_trigger, model, arguments = arguments);
    
    # add events to the model as reactions (handle roles)
    handle_events( events, model, arguments = arguments);
    
    ##### CLEANUP
    
    # add unprot ids
    if arguments.use_uniprot:
        cleanup_add_uniprot_information( model, entity_trigger, arguments = arguments);

    # remove reactions w/o reactants, products or modifiers
    if not arguments.remove_unconnected_reactions:
        cleanup_remove_unconnected_reactions( model, arguments = arguments);
    
    # remove unwanted species (the ones not involved in any reaction)
    if not arguments.remove_unconnected_species:
        cleanup_remove_unconnected_species( model, arguments = arguments);

    # complete reactions with reactants and products if none exist
    if arguments.complete_reactions:
        cleanup_complete_reactions( model, events, arguments = arguments);

    # remove empty compartments
    if not arguments.remove_empty_compartments:
        cleanup_remove_empty_compartments( model, arguments = arguments);
    
    return document

##################################################################
##### INITIALIZE MODEL, LOAD DATA etc
##################################################################

if __name__ == '__main__':
    
    
    CMD_ARGS = parse_arguments()
    
    logging.basicConfig( level = logging.INFO)
    logging.getLogger( "st2sbml")
    
    if CMD_ARGS.log_level == "ERROR":
        logging.getLogger( "st2sbml").setLevel( logging.ERROR);
    elif CMD_ARGS.log_level == "INFO":
        logging.getLogger( "st2sbml").setLevel( logging.INFO);
    elif CMD_ARGS.log_level == "WARNING":
        logging.getLogger( "st2sbml").setLevel( logging.WARNING);
    elif CMD_ARGS.log_level == "DEBUG":
        logging.getLogger( "st2sbml").setLevel( logging.DEBUG);
        
    ##### PARSE ANN or A1/A2
    if CMD_ARGS.path: # CMD_ARGS.path found
        if CMD_ARGS.use_ann:
            ANN_FILE_PATH = CMD_ARGS.path + CMD_ARGS.ann
            logging.getLogger( "st2sbml").info( "Processing %s", ANN_FILE_PATH)
            if not os.path.isfile( ANN_FILE_PATH):
                logging.getLogger( "st2sbml").error( "input {0} does not exist or is not a file".format( ANN_FILE_PATH))
                exit(1)
            TRIGGER, ENTITY_TRIGGER, EVENTS = parse_standoff.parse_ann( ANN_FILE_PATH);
        else:
            A1_FILE_PATH = CMD_ARGS.path + CMD_ARGS.a1
            A2_FILE_PATH = CMD_ARGS.path + CMD_ARGS.a2
            logging.getLogger( "st2sbml").info( "Processing %s a1/a2", CMD_ARGS.path)
            if not os.path.isfile( A1_FILE_PATH):
                logging.getLogger( "st2sbml").error( "Input %s does not exist or is not a file", A1_FILE_PATH)
                exit(1)
                
            if not os.path.isfile( A2_FILE_PATH):
                logging.getLogger( "st2sbml").error( "Input %s does not exist or is not a file", A2_FILE_PATH)
                exit(1)
            # parse the a1/a2 files
            TRIGGER, ENTITY_TRIGGER, EVENTS = parse_standoff.parse_a1_a2( A1_FILE_PATH, A2_FILE_PATH);
    else:
        ANN_FILE_PATH = CMD_ARGS.file
        logging.getLogger( "st2sbml").info( "Processing %s", ANN_FILE_PATH)
        if not os.path.isfile( ANN_FILE_PATH):
            logging.getLogger( "st2sbml").error( "input {0} does not exist or is not a file".format( ANN_FILE_PATH))
            exit(1)
        TRIGGER, ENTITY_TRIGGER, EVENTS = parse_standoff.parse_ann( ANN_FILE_PATH);
        
    DOCUMENT = create_document( trigger = TRIGGER, entity_trigger = ENTITY_TRIGGER, events = EVENTS, arguments = CMD_ARGS)
    
    ##### WRITE SBML FILE
    libsbml.writeSBMLToFile( DOCUMENT, CMD_ARGS.output)