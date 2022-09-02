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

Copyright 2021 RegulonDB

Permission to use, copy, modify, and/or distribute this software for any purpose with or without fee is hereby granted, provided that the above copyright notice and this permission notice appear in all copies.

THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
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


