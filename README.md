## TO-DO
- Issues with inconsistencies generally:
     - Detect number of fields in .res header to smartly grab results, e.g. how to figure out which data is missing in .res file
     - Consistently pull out correct data from directory name, e.g.               
```
#!bash
Ge+Te/190-0.07-1.75-1.75-PBE-00PBE 
Li+S/LiS-500-0.05-1.75-1.75-PBE-GO
```
e.g. Ge+Te/190-0.07-1.75-1.75-PBE-00PBE 
- Authentication with the database
- Provenance of structures (try to scrape CRSID from ssh)
- Automated tags conveying reliability and manual custom tags (e.g. 'me388 - miniproj 1')
- DB query functionality which mimics cryan at some level; fryan.py
- convert .castep file xc_functional to param acronym