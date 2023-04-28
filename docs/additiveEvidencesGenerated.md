---
title: "additiveEvidencesGenerated"
output:
  html_document: default
  pdf_document: default
date: "2023-04-28"
---


## Casos de Prueba de generación de additiveEvidences en datamart de Regulon

Este documento tiene como objetivo validar la generación de *additiveEvidences* en el datamart de Regulon. Dada su complejidad ha presentado algunos cambios y se han detectado algunas inconsistencias, por lo que estos casos de prueba permitiran verificar que la descarga es correcta.

Estos casos de prueba deben actualizarse en cada release, porque puede ocurrir que se agreguen nuevas evidencias, y el calculo de las additive evidences y confidence level puede ajustarse.

Ecocyc release: 26.5
ReguloDB pre-release 12.0-beta

-----

### VALIDACION de DESCARGA DE EVIDENCIAS

Casos recopilados  por: Andrés López Almazo  
Información de  identificadores de Ecocyc: Felipe Betancourt Figueroa 
Verificación de los casos: 
Validación:   

ESTADO:  

-----

Columnas de los registros mostrados en los casos de prueba para Regulon

```
# Columns:
# (1) tfId. TF id
# (2) tfName. TF name
# (3) cnfName. TF active conformation name
# (4) riId. Regulatory interaction (RI) identifier assigned by RegulonDB
# (5) riType. Regulatory interaction type [tf-promoter, tf-tu, tf-gene]
# (6) tfbsID. TF binding site (TFBS) identifier assigned by RegulonDB
# (7) tfbsLeft. TFBS left end position in the genome
# (8) tfbsright. TFBS right end position in the genome
# (9) strand. DNA strand where the TFBS is located
# (10) tfbsSeq. TFbS sequence (upper case)
# (11) function. Gene expression effect caused by the TF bound to the TFBS
# (12) promoterID. Promoter Identifier assigned by RegulonDB
# (13) promoterName. Promoter name
# (14) TSS. Transcription start site (+1) position in the genome
# (15) tfbsDistToPm. Relative distance from the center position of TFBS to the Transcription Start Site
# (16) firstGene. first transcribed gene name
# (17) tfrsDistTo1Gene. Relative distance from center position of TFBS to the start of first gene
# (18) target[tu|gene]. Transcription unit or gene (id:name) regulated by the TF
# (19) confidenceLevel. RI confidence level (Values: Confirmed, Strong, Weak)
# (20) tfrsEvidence. Evidence that supports the existence of the TFRS [EvidenceCode|EvidenceType(C:confirmed S:strong W:weak)]]
# (21) riEvidence. Evidence that supports the RI function [EvidenceCode|EvidenceType(C:confirmed S:strong W:weak)]
# (22) addEvidence. Additive Evidence [CV(EvidenceCode1/EvidenceCodeN)|Confidence Level]

```

| Status | 1)tfid| 2)tfName | 3)cnfName  | 4)riId RegulonDB| 5)riType    | 6)tfbsID         | 7)tfrsLeft | 8)tfrsright | 9)strand | 10)tfrsSeq | 11)function | 12)promoterID    | 13)promoterName | 14)TSS  | 15)tfbsDistToPm | 16)firstGene | 17)tfbsDistTo1Gene | 18)target[tu:gene] | 19)confidenceLevel | 20)tfbsEvidence                                                                                         | 21)riEvidence                                           | 22)addEvidence                                                                                                                                       |
|:---|:---|:---|:-----------|:---|:------------|:-----------------|:---|:---|:---------|:---|:---|:-----------------|:----------------|:--------|:----------------|:-------------|:-------------------|:-------------------|:-------------------|:--------------------------------------------------------------------------------------------------------|:--------------------------------------------------------|:-----------------------------------------------------------------------------------------------------------------------------------------------------|
|        | RDBECOLITFC00206 | PuuR     | PuuR       | RDBECOLIRIC00902 | tf-promoter | RDBECOLIBSC02314 | 1361085     |  1361085 | forward  | aagcgcagtaATGGCGATAATTTAGTCCACtttgtgagat | repressor | RDBECOLIPMC00069 | puuDp           | 1361052 | -      | -               | -25.5              | - | C                  | [IC:W], [EXP-IMP-SITE-MUTATION:S], [EXP-IDA-BINDING-OF-PURIFIED-PROTEINS:S]                             | [IC:W], [EXP-IEP-GENE-EXPRESSION-ANALYSIS:W]            | [AE(EXP-IMP-SITE-MUTATION/EXP-IDA-BINDING-OF-PURIFIED-PROTEINS):C]                                                                                   |
|        | RDBECOLITFC00006 | AtoC     | AtoC-P     | RDBECOLIRIC04777 | tf-promoter | RDBECOLIBSC03557 | 2323264     |  2323283 | forward  | ccaaaacttgCTATGCAGAAATTTGCACAGtgcgcaattt | activator | RDBECOLIPMC03436 | atoDp           | 2323412 | -      | -               | -173.5             | - | C                  | [COMP-HINF-SIMILAR-TO-CONSENSUS:W], [EXP-IMP-SITE-MUTATION:S], [EXP-IDA-BINDING-OF-PURIFIED-PROTEINS:S] | [EXP-IEP-GENE-EXPRESSION-ANALYSIS:W]                    | [AE(EXP-IMP-SITE-MUTATION/EXP-IDA-BINDING-OF-PURIFIED-PROTEINS):C]                                                                                   |
|        | RDBECOLITFC00111 | NorR     | NorR       | RDBECOLIRIC05142 | tf-promoter | RDBECOLIBSC00906 | 2832301     |  2832313 | forward  | tttgcctcacTGTCAATTTGACTatagatattg | activator | RDBECOLIPMC02745 | norVp           | 2832439 | -      | -               | -169               | - | C                  | [COMP-AINF-SIMILAR-TO-CONSENSUS:W], [EXP-GSELEX:W], [EXP-IDA-BINDING-OF-PURIFIED-PROTEINS:S]            | [EXP-IEP-GENE-EXPRESSION-ANALYSIS:W]                    | [AE(COMP-AINF-SIMILAR-TO-CONSENSUS/EXP-GSELEX/EXP-IEP-GENE-EXPRESSION-ANALYSIS):S],[AE(COMP-AINF-SIMILAR-TO-CONSENSUS/EXP-GSELEX/EXP-IEP-GENE-EXPRESSION-ANALYSIS/EXP-IDA-BINDING-OF-PURIFIED-PROTEINS):C]                                                              |
|        | RDBECOLITFC00108 | MetR     | MetR       | RDBECOLIRIC03430 | tf-promoter | RDBECOLIBSC03460 | 4223719     |  4223733 | forward  | tgagttaatgTTGAACAAATCTCATgttgcgtggt | activator | RDBECOLIPMC03363 | metHp1          | 4223783 | -      | -               | -102               | - | W                  | [COMP-HINF-SIMILAR-TO-CONSENSUS:W]                                                                      | [EXP-IEP-GENE-EXPRESSION-ANALYSIS:W]                    | -                                                                                                                                                    |
|        | RDBECOLITFC00135 | Nac     | Nac        | RDBECOLIRIC01569 | tf-tu       | RDBECOLIBSC02647 | 1404611     |  1404627 | -        | ttcacgtagcGATAGTTTTTACTTATCactaactgat | repressor | -                | -               | -       | -      | -               | -                  | - | S                  | [COMP-AINF-PATTERN-DISCOVERY:W],[EXP-CHIP-SEQ:W]                                                        | [EXP-IEP-RNA-SEQ:W]                    | [AE(COMP-AINF-PATTERN-DISCOVERY/EXP-CHIP-SEQ/EXP-IEP-RNA-SEQ):S]                                                                                     |