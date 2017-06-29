# Standoff-conversion
Python scripts for converting NLP standoff files to SBML and BioPAX.

Detailed conversion scheme at : http://sbnlp.github.io/standoff-conversion/

Read the complete paper here:


[Extracting Biological Pathway Models From NLP Event Representations](http://aclweb.org/anthology/W/W15/W15-3805.pdf) - Proceedings of the ACL 2015 Workshop on Biomedical Natural Language Processing (BioNLP’15)

#Prerequisites
* paxtools-4.2.1.jar; 
Available for download from http://sourceforge.net/projects/biopax/files/paxtools/

* Java to python integrator (JPype); See https://pypi.python.org/pypi/JPype1

* libSBML with python bindings; See http://sbml.org/Software/libSBML/docs/python-api/

#Usage
To  understand the different options available for the conversion from standoff to BioPAX/SBML

``` python st2sbml.py --help ```

```python st2biopax.py --help ```


# Copyright

Copyright [2015-2017] [Michael Spranger, Sucheendra K. Palaniappan, Samik Ghosh]

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this software except in compliance with the License.

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
