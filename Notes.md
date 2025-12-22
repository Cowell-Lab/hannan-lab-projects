# Kidney TCR project notes
project uuid: 4661bdd7-7498-412d-b218-f42044550b96

## VDJPipe pre-processing
The analyses document IDs are:


* `library 33 analysis document id`: 37afeb5b-ae5c-49a1-8c75-366f978c996a
* `library 34 analysis document id`: 91ca4180-37bd-4d93-aa7b-4425deee0b95

* `kidney_library_33_samples.json`: library 33 data, job id: 37afeb5b-ae5c-49a1-8c75-366f978c996a/9251db8a-f7ec-45a6-8f3b-79f25232154c-007
* `kidney_library_34_samples.json`: library 34 data, job id: 91ca4180-37bd-4d93-aa7b-4425deee0b95/6c3d5355-4429-4431-be5b-51162af8b849-007


## IgBlast annotation
Ran IGBlast on library 33 and 34 from VDJPipe output.

* `job-igblast-KidneyTCR.json`: library 33 and 34 data, job id: 0c31b9fc-9732-4eb0-acf9-67b015a0e12d-007
* `job-igblast-KidneyTCR.json`: library adaptive, job id: edfebb12-dd91-4a75-9827-338f6fee9f5c/1361103b-f356-4833-93ae-2a50cfd608c1-007
 
 
## RepCalc analysis
Ran Repcalc Job with the output from the first run of IGBlast. Ran this job with repcalc version 0.2

* `job-repcalc-KidneyTCR.json`: library 33 and 34 data, job id: 7ee927d9-c1b3-4fa2-a5bc-20bde3e773cd-007
  
Ran repcalc again with both igblast jobs. This is with repcalc version 0.5
* `job-repcalc-KidneyTCR.json`: library 33 and 34 data and adaptive samples, job id: 88a22c28-d155-4b97-8078-9ed86ec27e43-007

### RepCalc gene usage

TBD

### RepCalc gene combo usage

TBD

### RepCalc CDR3 length distribution

TBD

### RepCalc CDR3 sharing

TBD
