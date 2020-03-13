# significant_events/

This folder has been archived and can be accessed here:

The workflow involves either generating new data - that will be pulled into the relative directory

`sbas/data`

There is test data available here:

`https://github.com/adeslatt/sbas_test/releases/tag/rmats_final.gencode.v30`

To obtain that data do the following:

```
cd sbas/data
wget https://github.com/adeslatt/sbas_test/releases/download/rmats_final.gencode.v30/fromGTF.A3SS.txt
wget https://github.com/adeslatt/sbas_test/releases/download/rmats_final.gencode.v30/fromGTF.A5SS.txt
wget https://github.com/adeslatt/sbas_test/releases/download/rmats_final.gencode.v30/fromGTF.MXE.txt
wget https://github.com/adeslatt/sbas_test/releases/download/rmats_final.gencode.v30/fromGTF.RI.txt
wget https://github.com/adeslatt/sbas_test/releases/download/rmats_final.gencode.v30/fromGTF.SE.txt
```

Additionally, the significant result data from differential analysis may be obtained and pulled into
the relative directory:

`sbas/data/significant_events`

There is test data available here:

`https://github.com/adeslatt/sbas_test/releases/tag/GTExV6SignificantASTissueEvents.v1`

To obtain that data, do the following:

```
cd sbas/data
wget https://github.com/adeslatt/sbas_test/releases/download/GTExV6SignificantASTissueEvents.v1/significant_events.tar
mkdir significant_events
tar xvf significant_events
```
