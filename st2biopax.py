# Script for mapping standoff formats to biopax (see ... for more details)
# It combines information from
#  [1] T. Ohta, S. Pyysalo, J. Tsujii. From Pathways to Biomolecular Events: Opportunities and Challenges In Proceedings of the 2011 Workshop on Biomedical Natural Language Processing, ACL-HLT 2011, pages 105-113, Portland, Oregon, USA, June 23-24, 2011.
#      http://www.aclweb.org/anthology/W11-0214.pdf
#  [2] Pathway Curation (PC) task/BioNLP-ST 2013 http://2013.bionlp-st.org/tasks/pathway-curation
#  [3] BioPax level 3 documentation http://www.biopax.org/release/biopax-level3-documentation.pdf
#  [4] CellDesigner BioPax conversion tool

# to run all tests
# ls tests/*.ann | xargs -I {} python2.7 st2biopax.py --file {}  --output {}.output-biopax.owl --complete-interactions

import argparse
import copy
import inspect
import jpype
import logging
import os
import re
import sets

import parse_standoff
import uniprot_tools

ENTITY_CLASSES = { "cellular_component" : "org.biopax.paxtools.model.level3.PhysicalEntity",
                  "complex" : "org.biopax.paxtools.model.level3.Complex",
                  "dna" : "org.biopax.paxtools.model.level3.Dna",
                  "drug" : "org.biopax.paxtools.model.level3.PhysicalEntity",
                  "entity" : "org.biopax.paxtools.model.level3.PhysicalEntity",
                  "gene_or_gene_product" : "org.biopax.paxtools.model.level3.PhysicalEntity",
                  "gene_product" : "org.biopax.paxtools.model.level3.PhysicalEntity",
                  "gene" : "org.biopax.paxtools.model.level3.Gene",
                  "ion" : "org.biopax.paxtools.model.level3.PhysicalEntity",
                  "protein" : "org.biopax.paxtools.model.level3.Protein",
                  "receptor" : "org.biopax.paxtools.model.level3.PhysicalEntity",
                  "rna" : "org.biopax.paxtools.model.level3.Rna",
                  "simple_molecule" : "org.biopax.paxtools.model.level3.SmallMolecule",
                  "simple_chemical" : "org.biopax.paxtools.model.level3.SmallMolecule",
                  "tag" : "org.biopax.paxtools.model.level3.PhysicalEntity"};

# events: standoff -> biopax interaction
EVENT_CLASSES = {
    "acetylation" : { "type" : "org.biopax.paxtools.model.level3.BiochemicalReaction", 
                      "psi_mod_right" : ["acetylated residue", "MOD:00394"] },
    "activation" : { "type" : "org.biopax.paxtools.model.level3.Control" },
    "association" : { "type" : "org.biopax.paxtools.model.level3.ComplexAssembly" },
    "binding" : { "type" : "org.biopax.paxtools.model.level3.ComplexAssembly" },
    "catabolism" : { "type" : "org.biopax.paxtools.model.level3.Degradation"},
    "catalysis" : { "type" : "org.biopax.paxtools.model.level3.Catalysis" },
    "conversion" : { "type" : "org.biopax.paxtools.model.level3.Conversion" },
    "deacetylation" : { "type" : "org.biopax.paxtools.model.level3.BiochemicalReaction",
                       "psi_mod_left" : ["acetylated residue", "MOD:00394"]},
    "inactivation" : { "type" : "org.biopax.paxtools.model.level3.Control" },
    "degradation" : { "type" : "org.biopax.paxtools.model.level3.Degradation" }, 
    "demethylation" : { "type" : "org.biopax.paxtools.model.level3.BiochemicalReaction",
                        "psi_mod_left" : ["methylated residue", "MOD:00427"]},
    "dephosphorylation" : { "type" : "org.biopax.paxtools.model.level3.BiochemicalReaction", # sub of Conversion [2], BiochemicalReaction [4]
                           "psi_mod_left" : ["phosphorylated residue","MOD:00696"]},
    "deubiquitination" : { "type" : "org.biopax.paxtools.model.level3.BiochemicalReaction", 
                          "psi_mod_left" : ["ubiquitinylation residue", "MOD:00492"]},
    "dissociation" : { "type" : "org.biopax.paxtools.model.level3.ComplexAssembly"},
    "gene_expression" : { "type" : "org.biopax.paxtools.model.level3.TemplateReaction"},
    "localization" : { "type" : "org.biopax.paxtools.model.level3.Transport"},
    "methylation" : { "type" : "org.biopax.paxtools.model.level3.BiochemicalReaction", 
                     "psi_mod_right" : ["methylated residue","MOD:00427"] },
    "negative_regulation" : { "type" : "org.biopax.paxtools.model.level3.Control" },
    "pathway" : { "type" : "org.biopax.paxtools.model.level3.Interaction" },
    "phosphorylation" : { "type" : "org.biopax.paxtools.model.level3.BiochemicalReaction", # sub of Conversion/Interaction [2], BiochemicalReaction [4]
                         "psi_mod_right" : ["phosphorylated residue", "MOD:00696"] },
    "positive_regulation" : { "type" : "org.biopax.paxtools.model.level3.Catalysis" },
    "protein_catabolism" : { "type" : "org.biopax.paxtools.model.level3.Degradation" },
    "regulation" : { "type" : "org.biopax.paxtools.model.level3.Control" },
    "transcription" : { "type" : "org.biopax.paxtools.model.level3.TemplateReaction" },
    "translation"  : { "type" : "org.biopax.paxtools.model.level3.TemplateReaction" },

    "transport" : { "type" : "org.biopax.paxtools.model.level3.Transport" },

    "ubiquitination" : { "type" : "org.biopax.paxtools.model.level3.BiochemicalReaction",
                        "psi_mod_right" : ["ubiquitinylation residue", "MOD:00492"]},
};

def parse_arguments( args = None):

    parser = argparse.ArgumentParser(description='Converts BioNLP shared task format to BioPax (requires paxtools-4.2.1.jar or higher).')

    parser.add_argument( '--file',
                    action = "store",
                    dest = "file",
                    default = "test.ann",
                    help = "the file to convert (default is test.ann)")
    
    parser.add_argument( '--path',
                        action = "store",
                        dest = "path",
                        default = None,
                        help="the path to ann or a1/a2 file (includes the prefix of the file)")
    
    parser.add_argument( '--a1',
                        action = "store",
                        dest = "a1",
                        default = ".a1",
                        help = "the file ending for a1 files (a1 by default)")
    
    parser.add_argument( '--a2',
                        action = "store",
                        dest = "a2",
                        default = ".a2",
                        help = "the file ending for a2 files (a2 by default)")
    
    parser.add_argument( '--ann',
                        action = "store",
                        dest = "ann",
                        default = ".ann",
                        help = "the file ending for ann files (ann by default)")
    
    parser.add_argument('--use-ann',
                        action = "store_true",
                        dest = "use_ann",
                        default = True, 
                        help = "whether to use ann or a1/a2")
    
    parser.add_argument( '--output',
                        action = "store",
                        dest = "output",
                        default = "output-biopax.owl",
                        help = "the path to output file (biopax.xml by default)")
    
    parser.add_argument( '--xml-base',
                        action = "store",
                        dest = "xml_base",
                        default = "http://example.com/",
                        help = "xml base for the biopax model (http://example.com by default)")
    
    parser.add_argument( '--remove-unconnected-physical-entities',
                        action = "store_true",
                        dest = "remove_unconnected_physical_entities",
                        default = False,
                        help = "whether to keep unconnected physical entities (i.e. entities not used in interactions)")

    parser.add_argument( '--complete-interactions',
                        action = "store_true",
                        dest = "complete_interactions",
                        default = False,
                        help = "whether to complete interactions by adding participants.")

    parser.add_argument( '--remove-unconnected-interactions',
                        action = "store_true",
                        dest = "remove_unconnected_interactions",
                        default = False,
                        help = "whether to keep interactions without participants.")

    parser.add_argument( '--use-uniprot',
                        action = "store_true",
                        dest = "use_uniprot",
                        default = False,
                        help = "add uniprot information")
    
    parser.add_argument( '--uniprot-pickle',
                        action = "store",
                        dest = "uniprot_pickle",
                        default = None, 
                        help = "whether and where to get uniprot information from a pickle file (e.g. you can use\"uniprot.pickle\" which is generated by uniprot_tools.py)")

    parser.add_argument( '--log-level',
                    action = "store",
                    dest = "log_level",
                    default = "INFO",
                    help = "the log level for messages. Can be INFO, WARNING, DEBUG, ERROR (INFO by default).")
        
    return parser.parse_args( args = args)

DEFAULT_ARGUMENTS = parse_arguments()

def get_path( arguments):
    if arguments.path != None:
        return arguments.path
    else:
        return arguments.file
    
def get_java_class( string):
    return jpype.java.lang.Class.forName( string, True, jpype.java.lang.ClassLoader.getSystemClassLoader())

def get_compartment_rdf_id( compartment_id, arguments = DEFAULT_ARGUMENTS):
    return arguments.xml_base + "cellularLocationVocabulary_" + compartment_id

def create_compartment( compartment_id, compartment_name, model, arguments = DEFAULT_ARGUMENTS):
    compartment_rdf_id = get_compartment_rdf_id( compartment_id, arguments)
    if not model.getByID( compartment_rdf_id):
        cell_loc = model.addNew( get_java_class( "org.biopax.paxtools.model.level3.CellularLocationVocabulary"),
                                    compartment_rdf_id)
        cell_loc.addTerm( compartment_name)
        
def set_compartment( entity_rdf_id, compartment_id, model, arguments = DEFAULT_ARGUMENTS):
    physical_entity = model.getByID( entity_rdf_id)
    compartment_rdf_id = get_compartment_rdf_id( compartment_id, arguments)
    compartment = model.getByID( compartment_rdf_id);
    if compartment is None:
        logging.getLogger( "st2biopax").error( "{0} entity {1} compartment {2} not set, because not in the model".format( get_path( arguments), entity_rdf_id, compartment_id))
        return;
    elif physical_entity is None:
        logging.getLogger( "st2biopax").error( "{0} entity {1} compartment not set because entity not in the model".format( get_path( arguments), entity_rdf_id))
        return
    else:
        physical_entity.setCellularLocation( compartment);
        
def add_physical_entity( type, id, name, model):
    bio_pax_class_str = type;
    entity = model.addNew( get_java_class( type), id)
    entity.setStandardName( name);
    entity.addComment( "st2biopax: created " + entity.getRDFId())
    return entity

def add_physical_entity_for_trigger( trigger, model, entity_id = None, compartment = None, arguments = DEFAULT_ARGUMENTS):
    if entity_id is None:
        entity_id = arguments.xml_base + trigger.id
    entity = add_physical_entity( type = ENTITY_CLASSES[ trigger.type_lower],
                                id = entity_id,
                                name = trigger.text,
                                model = model)
    if not compartment is None:
        set_compartment( entity_id, compartment, model, arguments)
    return entity

def create_sequence_modification_vocabulary( id, term, xref_str, model, arguments = DEFAULT_ARGUMENTS):
    """ creates a sequence modification vocabulary or get it"""
    if model.getByID( id):
        return model.getByID( id)
    else:
        seq_mod_vocab = model.addNew( get_java_class( "org.biopax.paxtools.model.level3.SequenceModificationVocabulary"), id);
        seq_mod_vocab.addTerm( term)
        
        if not xref_str is None:
            xref_id = arguments.xml_base + "SequenceModificationVocabulary_" + xref_str
            xref = model.getByID( xref_id)
            if xref is None:
                xref = model.addNew( get_java_class( "org.biopax.paxtools.model.level3.UnificationXref"), xref_id)
            seq_mod_vocab.addXref( xref)
        return seq_mod_vocab

def get_sequence_modification_vocabulary( event_type, left, model, arguments = DEFAULT_ARGUMENTS):
    """ Retrieves or creates a sequence modification vocabulary """
    if left:
        event_mapping_name = "psi_mod_left"
        seq_mod_vocab_id = arguments.xml_base + "SequenceModificationVocabulary_" + event_type + "_left"
    else:
        event_mapping_name = "psi_mod_right"
        seq_mod_vocab_id = arguments.xml_base + "SequenceModificationVocabulary_" + event_type + "_right"
    
    seq_mod_vocab = model.getByID( seq_mod_vocab_id)    
    event_mapping = EVENT_CLASSES.get( event_type)
    if event_mapping:
        mapping = event_mapping.get( event_mapping_name)
    
    if seq_mod_vocab is None and event_mapping and mapping:
        seq_mod_vocab = create_sequence_modification_vocabulary( id = seq_mod_vocab_id,
                                                                term = mapping[0],
                                                                xref_str = mapping[1],
                                                                model = model,
                                                                arguments = arguments)
    return seq_mod_vocab

def get_modification_feature( id, model, arguments = DEFAULT_ARGUMENTS):
    """ Creates a modification feature for event or retrieves it """
    modification_feature = model.getByID( id + "_ModificationFeature")
    if modification_feature is None:
        modification_feature = model.addNew( get_java_class( "org.biopax.paxtools.model.level3.ModificationFeature"),
                                            id + "_ModificationFeature")
    return modification_feature;

def set_sequence_modification_feature( entity, interaction, left, model, arguments = DEFAULT_ARGUMENTS):
    vocab = get_sequence_modification_vocabulary( str( interaction.getStandardName()).lower(),
                                                 left = left,
                                                 model = model,
                                                 arguments = arguments)
    if vocab:
        modification_feature_id = interaction.getRDFId()
        if left:
            modification_feature_id += "_Left"
        else:
            modification_feature_id += "_Right"
        modification_feature = get_modification_feature( modification_feature_id, model, arguments = arguments);
        modification_feature.setModificationType( vocab)
        entity.addFeature( modification_feature)

def get_right( interaction, model, arguments = DEFAULT_ARGUMENTS, create_right_if_none_exist = True):
    """ Retrieves all right entities or creates one """
    # handle conversions
    if isinstance( interaction, get_java_class( "org.biopax.paxtools.model.level3.Conversion")):
        if interaction.getRight().size() > 0:
            return [ p for p in interaction.getRight().toArray()]
        elif create_right_if_none_exist:
            right_id = interaction.getRDFId() + "_Right_0"
            if interaction.getLeft().size() > 0 \
                and interaction.getLeft().toArray()[0].getStandardName():
                right_name = interaction.getLeft().toArray()[0].getStandardName()
            else:
                right_name = "Right"
            right = add_physical_entity( id = right_id,
                                        type = "org.biopax.paxtools.model.level3.PhysicalEntity",
                                        name = right_name,
                                        model = model);            
            interaction.addRight( right)
    
            set_sequence_modification_feature( right,
                                              interaction = interaction,
                                              left = False,
                                              model = model,
                                              arguments = arguments)
        return [ right]
    # handle template reactions
    elif isinstance( interaction, get_java_class( "org.biopax.paxtools.model.level3.TemplateReaction")):
        if interaction.getProduct().size() > 0:
            return [ p for p in interaction.getProduct().toArray()]
        else:
            product_id = interaction.getRDFId() + "_Product_0"
            if interaction.getTemplate() and interaction.getTemplate().getStandardName() is not None:
                product_name = interaction.getTemplate().getStandardName()
            else:
                product_name = "Product"
            product = add_physical_entity( id = product_id,
                                        type = "org.biopax.paxtools.model.level3.PhysicalEntity",
                                        name = product_name,
                                        model = model)
            interaction.addProduct( product)
    
        return [ product]
    else:
        return []
    
def get_left( interaction, model, arguments = DEFAULT_ARGUMENTS, create_left_if_none_exist = True):
    """ Retrieves all left entities or creates one """
    if interaction.getLeft().size() > 0:
        return [ p for p in interaction.getLeft().toArray()]
    elif interaction.getLeft().size() == 0 and create_left_if_none_exist:
        left_id = interaction.getRDFId() + "Left_0" + str( interaction.getLeft().size())
        left = add_physical_entity( id = left_id,
                                   type = "org.biopax.paxtools.model.level3.PhysicalEntity",
                                   name = "Left",
                                   model = model)
        interaction.addLeft( left);
        
        set_sequence_modification_feature( left,
                                          interaction = interaction,
                                          left = True,
                                          model = model,
                                          arguments = arguments)
        
        return [ left]
        # handle template reactions
    elif isinstance( interaction, get_java_class( "org.biopax.paxtools.model.level3.TemplateReaction")):
        if not interaction.getTemplate() is None:
            return [ interaction.getTemplate()]
        else:
            template_id = interaction.getRDFId() + "Template"
            template = add_physical_entity( id = template_id,
                                           type = "org.biopax.paxtools.model.level3.NucleicAcid",
                                           name = "Template",
                                           model = model);            
            interaction.setTemplate( template)

        return [ template]
    else:
        return []

def add_interaction( type, id, standard_name, name, model):
    interaction = model.addNew( get_java_class( type), id)
    interaction.addComment( "st2biopax: created " + interaction.getRDFId())
    if name:
        interaction.addName( name);
    interaction.setStandardName( standard_name);
    
    if isinstance( interaction, get_java_class( "org.biopax.paxtools.model.level3.Conversion")):
        direction_class = get_java_class( "org.biopax.paxtools.model.level3.ConversionDirectionType")
        interaction.setConversionDirection( direction_class.LEFT_TO_RIGHT);
    return interaction

def add_interaction_for_event( event, model, arguments = DEFAULT_ARGUMENTS):
    if EVENT_CLASSES.get( event.type_lower):
        return add_interaction( type = EVENT_CLASSES[event.type_lower]["type"],
                                id = arguments.xml_base + event.id,
                                standard_name = event.type,
                                name = event.trigger.text,
                                model = model)
    else:
        logging.getLogger( "st2biopax").error( "{0} event {1} not added to the model, because unknown".format( get_path( arguments), event))
        return None
    
##################################################################
##### HANDLE ENTITIEs
##################################################################
  
def handle_entities( entities, model, arguments = DEFAULT_ARGUMENTS):
    """ goes through a set of entities and adds them to the model as physical entities """
    for trigger in entities:
        if ENTITY_CLASSES.get( trigger.type_lower):
            add_physical_entity_for_trigger( trigger, model);
        else:
            logging.getLogger( "st2biopax").error( "{0} entity {1} not added to model, because no mapping to BioPAX".format( get_path( arguments), trigger))

##################################################################
##### HANDLE EVENTS
##################################################################
   
def handle_localization( event, model, arguments = DEFAULT_ARGUMENTS):
    to_loc = None;
    from_loc = None;
    at_loc = None;

    if len( event.get_roles("toloc")) == 1:
        to_loc = event.get_roles("toloc")[0]
    elif len( event.get_roles("toloc")) != 0:
        logging.getLogger( "st2biopax").warning("{0} event {1} ToLoc not handled because more than 1 ToLoc".format( get_path( arguments), event)) 
        if not isinstance( to_loc, parse_standoff.EntityTrigger):
            logging.getLogger( "st2biopax").warning("{0} event {1} ToLoc {2} not handled because it is not an entity".format( get_path( arguments), event, to_loc))
            to_loc = None;
    if len( event.get_roles("fromloc")) == 1:
        from_loc = event.get_roles("fromloc")[0]
        if not isinstance( from_loc, parse_standoff.EntityTrigger):
            logging.getLogger( "st2biopax").warning("{0} event {1} FromLoc {2} not handled because it is not an entity".format( get_path( arguments), event, from_loc))
            from_loc = None;
    elif len( event.get_roles("fromloc")) != 0:
        logging.getLogger( "st2biopax").warning("{0} event {1} FromLoc not handled because more than 1 FromLoc".format( get_path( arguments), event))
    if len( event.get_roles("atloc")) == 1:
        at_loc = event.get_roles("atloc")[0]
        if not isinstance( at_loc, parse_standoff.EntityTrigger):
            logging.getLogger( "st2biopax").warning("{0} event {1} AtLoc {2} not handled because it is not an entity".format( get_path( arguments), event, at_loc))
            at_loc = None;
    elif len( event.get_roles("atloc")) != 0:
        logging.getLogger( "st2biopax").warning("{0} event {1} AtLoc not handled because more than 1 AtLoc".format( get_path( arguments), event))
        
    if from_loc and at_loc:
        logging.getLogger( "st2biopax").warning("{0} event {1} not handled because FromLoc and AtLoc given at the same time.".format( get_path( arguments), event))
    elif to_loc and at_loc:
        logging.getLogger( "st2biopax").warning("{0} event {1} not handled because ToLoc and AtLoc given at the same time.".format( get_path( arguments), event))
    elif from_loc is None and to_loc is None and at_loc is None:
        logging.getLogger( "st2biopax").warning("{0} event {1} not handled because no AtLoc, FromLoc or ToLoc given".format( get_path( arguments), event))
    else:
        if at_loc: # at_loc
            create_compartment( at_loc.id, at_loc.text, model = model, arguments = arguments)
            for theme in event.get_roles("Theme"):
                if theme.type_lower in parse_standoff.ENTITY_TRIGGER_TYPES:
                    set_compartment( entity_rdf_id = arguments.xml_base + theme.id, compartment_id = at_loc.id, model = model, arguments = arguments)
                else:
                    logging.getLogger( "st2biopax").warning("{0} event {1} theme {2} not handled because not a bp:PhysicalEntity but an bp:Interaction".format( get_path( arguments), event, theme))
        else: # from_loc and to_loc
            interaction = add_interaction_for_event( event, model, arguments = arguments);
    
            for theme in event.get_roles("Theme"):
                if not isinstance( theme, parse_standoff.EntityTrigger):
                    logging.getLogger( "st2biopax").warning("{0} event {1} theme {2} not handled because it is not an entity".format( get_path( arguments), event, theme))
                elif theme.type_lower in parse_standoff.ENTITY_TRIGGER_TYPES:
                    # add theme to interaction
                    theme_entity = model.getByID( arguments.xml_base + theme.id)
                    interaction.addLeft( theme_entity)
                    if from_loc:
                        create_compartment( from_loc.id, from_loc.text, model = model, arguments = arguments)
                        set_compartment( theme_entity.getRDFId(), compartment_id = from_loc.id, model = model, arguments = arguments);
                    if to_loc:
                        create_compartment( to_loc.id, to_loc.text, model = model, arguments = arguments)
                        to_entities = get_right( interaction, model = model, arguments = arguments)
                        for to_entity in to_entities:
                            set_compartment( to_entity.getRDFId(), compartment_id = to_loc.id, model = model)
                        
def handle_gene_expression( event, model, arguments = DEFAULT_ARGUMENTS):                 
    interaction = add_interaction_for_event( event, model, arguments = arguments);
    # for translation proteins are products (everything else is modifier)
    if event.type_lower == "translation":
        for theme in event.get_roles( "theme"):
            theme_entity = model.getByID( arguments.xml_base + theme.id)
            if isinstance( theme_entity, get_java_class( "org.biopax.paxtools.model.level3.NucleicAcid")):
                interaction.setTemplate( theme_entity);
            else:
                interaction.addProduct( theme_entity);
    # for gene_expression and transcription - rna and proteins are products
    else:    
        for theme in event.get_roles( "theme"):
            theme_entity = model.getByID( arguments.xml_base + theme.id)
            if isinstance( theme_entity, get_java_class( "org.biopax.paxtools.model.level3.Dna")):
                interaction.setTemplate( theme_entity);
            else:
                interaction.addProduct( theme_entity);

def complete_regulation_controlled( regulation, controlled, model, arguments = DEFAULT_ARGUMENTS):
    if regulation.getStandardName().lower() in [ "activation", "positive_regulation"]:
        id = arguments.xml_base + "SequenceModificationVocabulary_" + "Activated"
        seq_mod_vocab = create_sequence_modification_vocabulary( id = id,
                                                              term = "activated",
                                                              xref_str = None,
                                                              model = model)
        if controlled.getRight().size() == 0:
            get_right( controlled, model = model, arguments = arguments)
        for r in controlled.getRight().toArray():
            modification_feature = get_modification_feature( r.getRDFId(), model, arguments = arguments);
            modification_feature.setModificationType( seq_mod_vocab)
            r.addFeature( modification_feature)
            
    elif regulation.getStandardName().lower() in [ "inactivation", "negative_regulation"]:
        id = arguments.xml_base + "SequenceModificationVocabulary_" + "Inactivated"
        seq_mod_vocab = create_sequence_modification_vocabulary( id = id,
                                                              term = "inactivated",
                                                              xref_str = None,
                                                              model = model)
        if controlled.getRight().size() == 0:
            get_right( controlled, model = model, arguments = arguments)
        for r in controlled.getRight().toArray():
            modification_feature = get_modification_feature( r.getRDFId(), model, arguments = arguments);
            modification_feature.setModificationType( seq_mod_vocab)
            r.addFeature( modification_feature)

def handle_regulation( event, model, arguments = DEFAULT_ARGUMENTS):
    """ handles regulation events """
    unhandled_roles = []
    
    # we add interaction all themes and causes are entities
    if all( [ isinstance( theme, parse_standoff.EntityTrigger) for theme in event.get_roles( "theme")])\
        and all( [ isinstance( cause, parse_standoff.EntityTrigger) for cause in event.get_roles( "cause")]):
        interaction = add_interaction_for_event( event, model, arguments = arguments)
        control_type_class = get_java_class( "org.biopax.paxtools.model.level3.ControlType")
        
            
        # set control type
        if event.type_lower in ["positive_regulation", "activation"]:
            interaction.setControlType( control_type_class.ACTIVATION)
        elif event.type_lower in ["negative_regulation", "inactivation"]:
            interaction.setControlType( control_type_class.INHIBITION)

        for theme in event.get_roles("theme"):
            theme_id = arguments.xml_base + theme.id
            theme_entity = model.getByID( theme_id)

            proxy_interaction = add_interaction( type = "org.biopax.paxtools.model.level3.BiochemicalReaction",
                                                 id = arguments.xml_base + event.id + "_" + theme.id,
                                                 standard_name = "BiochemicalReaction",
                                                 name = None,
                                                 model = model)
            proxy_interaction.addLeft( theme_entity)
            interaction.addControlled( proxy_interaction);
            complete_regulation_controlled( interaction, proxy_interaction, model = model, arguments = arguments)

        for cause in event.get_roles( "cause"):
            interaction.addController( model.getByID( arguments.xml_base + cause.id))

    else:
        unhandled_roles = event.roles;

    # return unhandled roles
    if len(unhandled_roles) > 0:
        return { event.id : unhandled_roles }
    else:
        return None;

def resolve_cause_ids( cause, model, arguments = DEFAULT_ARGUMENTS):
    """ takes an entity or event (which is the cause of a regulation event and retrieves the actual cause id in RDF)
        in case there is none, it invents one """
    # cause can be an entity -> cause is the entity itself
    if isinstance( cause, parse_standoff.EntityTrigger):
        return [ arguments.xml_base + cause.id]
    # cause can be an event which is Regulation and not in model
    # -> cause is the cause of that event
    elif  isinstance( cause, parse_standoff.Event) \
        and cause.type_lower in [ "activation", "inactivation", "positive_regulation", "negative_regulation", "regulation", "catalysis"] \
        and len( cause.get_roles("cause")) > 0:
        # find the causes of the cause event
        results = []
        for c in cause.get_roles("cause"):
            results.extend( resolve_cause_ids( c, model, arguments = arguments))
        return results
    elif isinstance( cause, parse_standoff.Event) \
        and cause.type_lower in [ "activation", "inactivation", "positive_regulation", "negative_regulation", "regulation", "catalysis"]:
        cause_entity = add_physical_entity( id = cause.id + "_Cause_0",
                                           type = "org.biopax.paxtools.model.level3.PhysicalEntity",
                                           name = "Cause",
                                           model = model)
        return [ cause_entity.getRDFId()]
    elif  cause.type_lower in ["gene_expression", "transcription", "translation"]:
        interaction = model.getByID( arguments.xml_base + cause.id)
        return [ p.getRDFId() for p in interaction.getProduct()]
    # cause can be an event which is in the model
    # -> then the cause is actually the product of that event
    elif model.getByID( arguments.xml_base + cause.id):
        interaction = model.getByID( arguments.xml_base + cause.id)
        return [ p.getRDFId() for p in get_right( interaction, model)]
    # cannot handle other causes
    else:
        return [];

def resolve_themes( theme, model, arguments = DEFAULT_ARGUMENTS):
    """ takes an event (which is the theme of a regulation event and retrieves it's event) """
    # theme can be a Regulation event - this cannot be the theme, because it cannot be the controlled in BioPax
    if theme.type_lower in ["regulation", "negative_regulation", "positive_regulation", "catalysis"]:
        themes = []
        for theme in theme.get_roles( "theme"):
            themes.extend( resolve_themes( theme, model))
        return themes
    ## finds the theme can either be an entity or an event that exists in the model
    # cause can be an entity or an event
    elif  model.getByID( arguments.xml_base + theme.id):
        return [theme]
    else: 
        return [];

def handle_events( events, model, arguments = DEFAULT_ARGUMENTS):
    """ Processes all events in events and adds them to the model """
    unhandled_events = {}

    # pass 1 add events
    for event in events.values():
        if EVENT_CLASSES.get( event.type_lower) is None:
            # event is unknown
            logging.getLogger( "st2biopax").warning( "{0} event {1} unhandled, because unknown event".format( get_path( arguments), event))
        if event.type_lower == "pathway":
            pass # all entities are already added
        elif event.type_lower in [ "localization", "transport"]: 
            handle_localization( event, model, arguments = arguments);
        # handle regulation events (special handling)
        elif event.type_lower in [ "regulation", "positive_regulation", "negative_regulation", "activation", "inactivation", "catalysis"]:
            unhandled_event = handle_regulation( event, model, arguments = arguments);
            if not unhandled_event is None:
                unhandled_events.update( unhandled_event);
        elif event.type_lower in [ "gene_expression", "transcription", "translation"]:
            handle_gene_expression( event, model, arguments)
        # not all roles are entities
        elif not all( [ isinstance( role[1] , parse_standoff.EntityTrigger) for role in event.roles]):
            logging.getLogger( "st2biopax").warning( "{0} event {1} unhandled. Some roles are events, which is not allowed for this event type".format( get_path( arguments), event))
        # everything else: Conversion, Acetylation, Deacetylation, Demethylation, Dephosphorylation, Deubiquitination, Methylation, Phosphorylation, Ubiquitination
        else: 
            # add interaction
            interaction = add_interaction_for_event( event, model, arguments = arguments);
            #interaction.setSpontaneous(True)
            
            # handle products -> add as product
            for product in event.get_roles( "product"):
                right_physical_entity = model.getByID( arguments.xml_base + product.id)
                interaction.addRight( right_physical_entity)
                set_sequence_modification_feature( right_physical_entity,
                                                  interaction = interaction,
                                                  left = False,
                                                  model = model,
                                                  arguments = arguments)
            # handle themes -> add as reactants
            for theme in event.get_roles( "theme"):
                left_physical_entity = model.getByID( arguments.xml_base + theme.id);
                interaction.addLeft( left_physical_entity)
                set_sequence_modification_feature( left_physical_entity,
                                                  interaction = interaction,
                                                  left = True,
                                                  model = model,
                                                  arguments = arguments)
            # handle comp -> add as reactants
            for comp in event.get_roles( "complex"):
                interaction.addLeft( model.getByID( arguments.xml_base + comp.id))
            # handle Participant -> add as product
            for part in event.get_roles( "participant"):
                interaction.addLeft( model.getByID( arguments.xml_base + part.id))
            # # handle causes -> add as modifiers
            for cause in event.get_roles( "cause"):
                cause_biopax_element = model.getByID( arguments.xml_base + cause.id)
                catalysis_interaction = add_interaction( type = "org.biopax.paxtools.model.level3.Catalysis",
                                                         id = arguments.xml_base + event.id + "_Cause_" + cause.id + "_Catalysis",
                                                         standard_name = "Catalysis",
                                                         name = None,
                                                         model = model);
                catalysis_interaction.addController( cause_biopax_element)
                catalysis_interaction.addControlled( interaction)

            ## handle site -> add as modification feature
            for site in event.get_roles( "site"):
                # set FeatureLocation/setFeatureLocationType on the modification feature
                sequence_site = model.addNew( get_java_class( "org.biopax.paxtools.model.level3.SequenceSite"),
                                             arguments.xml_base + "SequenceSite_" + event.id + site.id)
                sequence_site.addComment( site.text)
                sequence_region_vocabulary = model.addNew( get_java_class( "org.biopax.paxtools.model.level3.SequenceRegionVocabulary"),
                                                          arguments.xml_base + "SequenceRegionVocabulary_" + event.id + "_Site_" + site.id)
                sequence_region_vocabulary.addTerm( site.text)
                left = get_sequence_modification_vocabulary( event.type_lower, left = True, model = model, arguments = arguments)
                if left:
                    modification_feature = get_modification_feature( interaction.getRDFId() + "_Left", model, arguments = arguments);
                else:
                    modification_feature = get_modification_feature( interaction.getRDFId() + "_Right", model, arguments = arguments);
                if modification_feature:
                    modification_feature.setFeatureLocation( sequence_site)
                    modification_feature.setFeatureLocationType( sequence_region_vocabulary)

            # set spontaneous
            #if event.type in [ "Binding", "Dissociation"]: 
            #    interaction.setSpontaneous(True)

            # check if there are any unhandled roles
            for unhandled_role in set([role[0] for role in event.roles]).difference(["theme", "product", "site", "participant", "complex", "cause"]):
                logging.getLogger( "st2biopax").warning( "{0} event {1} role {2} not handled (not allowed)".format( get_path( arguments), event, unhandled_role))


    # pass 2 add all unhandled events as interactions (only regulation events end up here)
    interactions = []
    for event_id in unhandled_events:
        event = events[event_id];
        interaction = add_interaction_for_event( event, model = model, arguments = arguments);
        control_type_class = get_java_class( "org.biopax.paxtools.model.level3.ControlType")
        if event.type_lower in ["positive_regulation", "activation"]:
            interaction.setControlType( control_type_class.ACTIVATION)
        elif event.type_lower in ["negative_regulation", "inactivation"]:
            interaction.setControlType( control_type_class.INHIBITION)
        interactions.append( interaction)

    # pass 3 add themes and causes for all unhandled events
    for idx, event_id in enumerate( unhandled_events):
        event = events[event_id]
        interaction = interactions[idx]
        
        
        for theme in event.get_roles( "theme"):
            if isinstance( theme, parse_standoff.EntityTrigger):
                proxy_interaction = add_interaction( type =  "org.biopax.paxtools.model.level3.BiochemicalReaction",
                                                     id = arguments.xml_base + event.id + "_" + theme.id,
                                                     standard_name = "BiochemicalReaction",
                                                     name = None,
                                                     model = model)
                theme_id = arguments.xml_base + theme.id
                theme_entity = model.getByID( theme_id)
                proxy_interaction.addLeft( theme_entity)
                interaction.addControlled( proxy_interaction);
                complete_regulation_controlled( interaction, proxy_interaction, model = model, arguments = arguments)
            else:
                real_themes = resolve_themes( theme, model)
                if real_themes is None:
                    logging.getLogger( "st2biopax").warning( "{0} event {1} theme {2} not handled, because could not be resolved".format( get_path( arguments), event, cause))
                else:
                    for real_theme in real_themes:
                        theme_biopax_element = model.getByID( arguments.xml_base + theme.id)
                        try:
                            interaction.addControlled( theme_biopax_element)
                        except:
                            logging.getLogger( "st2biopax").error( "{0} event {1} could not add controlled {2}".format( get_path( arguments), event, theme.id))
        for cause in event.get_roles( "Cause"):
            resolved_cause_ids = resolve_cause_ids( cause, model = model, arguments = arguments)
            for cause_id in resolved_cause_ids:
                cause_biopax_element = model.getByID( cause_id)
                interaction.addController( cause_biopax_element)
            if len( resolved_cause_ids) == 0:
                logging.getLogger( "st2biopax").warning( "{0} event {1} cause {2} not handled, because could not be resolved to an entity".format( get_path( arguments), event, cause))
        
        
                
##################################################################
##### POSTPROCESSING OPERATIONS
##################################################################

def cleanup_add_uniprot_information( model, entities, arguments = DEFAULT_ARGUMENTS):
    for entity in entities:
        physical_entity = model.getByID( arguments.xml_base + entity.id)
        if not physical_entity is None:
            uniprot_data = uniprot_tools.get_uniprot_data( entity.text, uniprot_pickle = arguments.uniprot_pickle) 
            if len( uniprot_data.keys()) == 0:
                logging.getLogger( "st2biopax").warning( "{0} no uniprot information found for entity {1}".format( get_path( arguments), entity))
            else:
                if not uniprot_data.get( "official_name") is None:
                    physical_entity.setStandardName( uniprot_data.get("official_name"))
                #if not uniprot_data.get("gene_id") is None:
                #    add_note( "Gene ID: {0}".format( uniprot_data.get("gene_id")), physical_entity, arguments = arguments);
                if not uniprot_data.get( "gene_name") is None:
                    physical_entity.addName( uniprot_data.get( "gene_name"));
                for alternate_name in uniprot_data.get( "alternate_names"):
                    physical_entity.addName( alternate_name);
                for gene_name_synonym in uniprot_data.get( "gene_name_synonyms"):
                    physical_entity.addName( gene_name_synonym);
                for uniprot_id in uniprot_data[ "uniprot_ids"]:
                    uniprot_id_str = "http://identifiers.org/uniprot/" + str( uniprot_id)
                    if model.getByID( uniprot_id_str) != None:
                        xref = model.getByID( uniprot_id_str)
                    else:
                        xref_class_str = "org.biopax.paxtools.model.level3.UnificationXref"
                        xref_class = get_java_class(  xref_class_str)
                        xref = model.addNew( xref_class, uniprot_id_str)
                        xref.setDb( "uniprot")
                    physical_entity.addXref( xref)

def cleanup_remove_unconnected_physical_entities( model, arguments = DEFAULT_ARGUMENTS):
    interactions = model.getObjects( get_java_class( "org.biopax.paxtools.model.level3.Interaction")).toArray()
    entity_type = get_java_class( "org.biopax.paxtools.model.level3.PhysicalEntity")
    entities = sets.Set( model.getObjects( entity_type).toArray())
    participants = sets.Set()
    for interaction in interactions:
        participants.update( [ participant for participant in interaction.getParticipant().toArray() if isinstance( participant, entity_type)]);
    
    non_participants = participants.difference( entities)
    
    for non_participant in non_participants:
        model.remove( non_participant)

def cleanup_remove_unconnected_interactions( model, arguments = DEFAULT_ARGUMENTS):
    # remove interactions w/o participants
    interactions = model.getObjects( get_java_class( "org.biopax.paxtools.model.level3.Interaction")).toArray()
    for interaction in interactions:
        if isinstance( interaction, get_java_class( "org.biopax.paxtools.model.level3.Degradation")):
            if interaction.getLeft().size() == 0:
                model.remove( interaction)
        elif isinstance( interaction, get_java_class( "org.biopax.paxtools.model.level3.Conversion")):
            if interaction.getLeft().size() == 0 and interaction.getRight().size() == 0:
                model.remove( interaction)
        elif isinstance( interaction, get_java_class( "org.biopax.paxtools.model.level3.TemplateReaction")):
            if interaction.getProduct().size() == 0 and interaction.getTemplate() == None:
                model.remove(interaction)
        elif isinstance( interaction, get_java_class( "org.biopax.paxtools.model.level3.Control")):
            if interaction.getControlled().size() and interaction.getController().size() == 0:
                model.remove( interaction)

def cleanup_complete_interactions( model, arguments = DEFAULT_ARGUMENTS):
    # complete interactions with reactants and products if none exist
    interactions = model.getObjects( get_java_class( "org.biopax.paxtools.model.level3.Interaction")).toArray()
    for interaction in interactions:
        if isinstance( interaction, get_java_class( "org.biopax.paxtools.model.level3.Degradation")):
            get_left( interaction, model = model, arguments = arguments);
        elif isinstance( interaction, get_java_class( "org.biopax.paxtools.model.level3.Conversion")):
            get_left( interaction, model = model, arguments = arguments)
            get_right( interaction, model = model, arguments = arguments)
        elif isinstance( interaction, get_java_class( "org.biopax.paxtools.model.level3.TemplateReaction")):
            if interaction.getTemplate() == None:
                template = add_physical_entity( type = "org.biopax.paxtools.model.level3.NucleicAcid",
                                                id = interaction.getRDFId() + "_template",
                                                name = "Template",
                                                model = model)
                interaction.setTemplate( template)
            if interaction.getProduct().size() == 0:
                if interaction.getTemplate() \
                 and interaction.getTemplate().getStandardName():
                    product_name = interaction.getTemplate().getStandardName()
                else:
                    product_name = "Product"
                product = add_physical_entity( type = "org.biopax.paxtools.model.level3.PhysicalEntity",
                                               id = interaction.getRDFId() + "_product",
                                               name =  product_name,
                                               model = model)
                interaction.addProduct( product);
        elif isinstance( interaction, get_java_class( "org.biopax.paxtools.model.level3.Control")):
            if interaction.getControlled().size() == 0:
                controlled = add_interaction( type = "org.biopax.paxtools.model.level3.BiochemicalReaction",
                                             id = interaction.getRDFId() + "_controlled",
                                             standard_name = "BiochemicalReaction",
                                             name = "Controlled",
                                             model = model);
                interaction.addControlled( controlled)
                get_left( controlled, model = model, arguments = arguments)
                get_right( controlled, model = model, arguments = arguments)
                complete_regulation_controlled( interaction, controlled, model = model, arguments = arguments)
            if interaction.getController().size() == 0:
                controller = add_physical_entity( type = "org.biopax.paxtools.model.level3.PhysicalEntity",
                                                  id = interaction.getRDFId() + "_controller",
                                                  name = "Controlled",
                                                  model = model)
                interaction.addController( controller)


##################################################################
##### CREATE MODEL
##################################################################

def initialize():
    """ call this only once """
    jpype.startJVM( jpype.getDefaultJVMPath(), "-ea", "-Xmx1g", "-Djava.class.path=paxtools-4.2.1.jar")
    
def create_model( entities, events, arguments = DEFAULT_ARGUMENTS):
    #call this to initialize use of Java
    
    #get the paxtools root package as a shortcut
    pax_pkg = jpype.JPackage("org.biopax.paxtools")

    #create a new BioPAX L3 factory
    l3factory = pax_pkg.model.BioPAXLevel.L3.getDefaultFactory()
    #create a new empty BioPAX model
    model = l3factory.createModel()
    #(highly recommended) use an xml base (URI prefix for elements we create)
    model.setXmlBase( arguments.xml_base)
    
    handle_entities( entities = entities, model = model, arguments = arguments)

    handle_events( events = events, model = model, arguments = arguments)
    
    ##### CLEANUP

    # add unprot ids    
    if arguments.use_uniprot:
        cleanup_add_uniprot_information( model, entities, arguments = arguments);

    # remove interactions w/o reactants, products or modifiers
    if not arguments.remove_unconnected_interactions:
        cleanup_remove_unconnected_interactions( model, arguments = arguments);
    
    # remove unwanted physical_entities (the ones not involved in any interaction)
    if not arguments.remove_unconnected_physical_entities:
        cleanup_remove_unconnected_physical_entities( model, arguments = arguments);

    # complete interactions with reactants and products if none exist
    if arguments.complete_interactions:
        cleanup_complete_interactions( model = model, arguments = DEFAULT_ARGUMENTS);

    # repair the model
    model.repair()
    return model

##################################################################
##### MAIN
##################################################################

if __name__ == '__main__':
    
    CMD_ARGS = parse_arguments()
    
    logging.basicConfig( level = logging.INFO)
    logging.getLogger( "st2biopax")

    # setup log levels    
    if CMD_ARGS.log_level == "ERROR":
        logging.getLogger( "st2biopax").setLevel( logging.ERROR);
    elif CMD_ARGS.log_level == "INFO":
        logging.getLogger( "st2biopax").setLevel( logging.INFO);
    elif CMD_ARGS.log_level == "WARNING":
        logging.getLogger( "st2biopax").setLevel( logging.WARNING);
    elif CMD_ARGS.log_level == "DEBUG":
        logging.getLogger( "st2biopax").setLevel( logging.DEBUG);
    
    # input, output files and global arguments
    OUTPUT_FILE = CMD_ARGS.output
    USE_UNIPROT = CMD_ARGS.use_uniprot # can True or False, if true uses uniprot_tools to add uniprot information for proteins
    USE_UNIPROT_PICKLE_FILE = CMD_ARGS.uniprot_pickle # use uniprot.pickle or online connection
    
    ##### PARSE ANN or A1/A2
    if CMD_ARGS.path: # CMD_ARGS.path found
        if CMD_ARGS.use_ann:
            ANN_FILE_PATH = CMD_ARGS.path + CMD_ARGS.ann
            logging.getLogger( "st2biopax").info( "Processing %s", ANN_FILE_PATH)
            if not os.path.isfile( ANN_FILE_PATH):
                logging.getLogger( "st2biopax").error( "input {0} does not exist or is not a file".format( ANN_FILE_PATH))
                exit(1)
            TRIGGER, ENTITY_TRIGGER, EVENTS = parse_standoff.parse_ann( ANN_FILE_PATH);
        else:
            A1_FILE_PATH = CMD_ARGS.path + CMD_ARGS.a1
            A2_FILE_PATH = CMD_ARGS.path + CMD_ARGS.a2
            logging.getLogger( "st2biopax").info( "Processing %s a1/a2", CMD_ARGS.path)
            if not os.path.isfile( A1_FILE_PATH):
                logging.getLogger( "st2biopax").error( "Input %s does not exist or is not a file", A1_FILE_PATH)
                exit(1)
                
            if not os.path.isfile( A2_FILE_PATH):
                logging.getLogger( "st2biopax").error( "Input %s does not exist or is not a file", A2_FILE_PATH)
                exit(1)
            # parse the a1/a2 files
            TRIGGER, ENTITY_TRIGGER, EVENTS = parse_standoff.parse_a1_a2( A1_FILE_PATH, A2_FILE_PATH);
    else:
        ANN_FILE_PATH = CMD_ARGS.file
        logging.getLogger( "st2biopax").info( "Processing %s", ANN_FILE_PATH)
        if not os.path.isfile( ANN_FILE_PATH):
            logging.getLogger( "st2biopax").error( "input {0} does not exist or is not a file".format( ANN_FILE_PATH))
            exit(1)
        TRIGGER, ENTITY_TRIGGER, EVENTS = parse_standoff.parse_ann( ANN_FILE_PATH);

    ############################
    # start jpype
    initialize();
    
    ############################
    # generating biopax model
    MODEL = create_model( entities = ENTITY_TRIGGER, events = EVENTS, arguments = CMD_ARGS)
    
    #export the model to a BioPAX OWL file
    PAX_PKG = jpype.JPackage("org.biopax.paxtools")
    JAVA_IO = jpype.JPackage("java.io")
    IO = PAX_PKG.io.SimpleIOHandler( PAX_PKG.model.BioPAXLevel.L3)
    OUTPUT_FILE = CMD_ARGS.output
    FILE_OS = JAVA_IO.FileOutputStream( OUTPUT_FILE)
    IO.convertToOWL( MODEL, FILE_OS)
    FILE_OS.close()
    
    #end use of jpype - docs say you can only do this once, so all java must be run before calling this
    #jpype.shutdownJVM() 
