# RegulonDB Datamarts

# Description

This software contains all needed to extract all data for RegulonDB 
Datamarts Collections. At this moment the collections that can be 
extracted are the next listed:
- Gene
- Regulon
- Operon
- Drawing Traces Tool (DTT)
- Gene Coexpression
- Sigmulon
- SRNA
- Regulatory Network

# Motivation
With the reengineering of RegulonDB passing from relational model to an documental model that contains all data of collection named now as datamart, with this app can be possible to get data from RegulonDB Multigenomic model and build the datamarts for the GraphQL API Query Services.


# System requirements

- Python ^3.9.7
- RegulonDB Multigenomic API Services
- 8GB RAM
- 50GB Space on Disk

# Install

Clone this repo with
```
git clone https://github.com/regulondbunam/RegulonDB-Datamarts
```

# Quick start

First install RegulonDB Multigenomic API from its Repo in this [link](https://github.com/regulondbunam/multigenomic-api)

Then you can use the following command to get all Datamarts Extraction files in JSON
```
python3 __main__.py
```

The extraction might take a while, the results will be in lib/data

# Project website 

[NOT DEFINED]

# License

Copyright 2023 RegulonDB

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

# Support contact information

[It should be clear where to go for support, for example a contact e-mail address]

# Software quality checklist

**Accessibility**

- [ ] Unique DOI [identifier](http://....) (Please update identifier and link)
- [ ] Version control system

**Documentation**

- [X] README file

**Learnability**

- [X] Quick start

**Buildability**

- [ ] INSTALL file

**Identity**

- [ ] Website

**Copyright & Licensing**

- [X] LICENSE file

**Portability**

- [ ] Multiple platforms
- [ ] Browsers

**Supportability**

- [ ] E-mail address
- [ ] Issue tracker
- [ ] Slack
- [ ] Gitter

**Analysability**

- [ ] Source code structured
- [ ] Sensible names
- [ ] Coding standards - [style guides](http://google.github.io/styleguide/)

**Changeability**

- [ ] CONTRIBUTING file
- [ ] Code of Conduct file
- [ ] Code changes, and their authorship, publicly visible

**Reusability**

- [ ] Source code set up in a modular fashion

**Security & Privacy**

- [ ] Passwords must never be stored in unhashed form


