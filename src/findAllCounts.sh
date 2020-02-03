#!/bin/bash
# gt 1000
awk -F ',' '{ for (i=1;i<=NF;i++) if ($i>1000) count[i]+=1} END{for (i in count) printf("%d ", count[i]); printf("\n") }' rmats.se.jc.sjc.male.txt > rmats.se.jc.sjc.male.rows.gt.1000.txt
awk -F ',' '{ for (i=1;i<=NF;i++) if ($i>1000) count[i]+=1} END{for (i in count) printf("%d ", count[i]); printf("\n") }' rmats.se.jc.sjc.female.txt > rmats.se.jc.sjc.female.rows.gt.1000.txt
awk -F ',' '{ for (i=1;i<=NF;i++) if ($i>1000) count[i]+=1} END{for (i in count) printf("%d ", count[i]); printf("\n") }' rmats.se.jc.ijc.female.txt > rmats.se.jc.ijc.female.rows.gt.1000.txt
awk -F ',' '{ for (i=1;i<=NF;i++) if ($i>1000) count[i]+=1} END{for (i in count) printf("%d ", count[i]); printf("\n") }' rmats.se.jc.ijc.male.txt > rmats.se.jc.ijc.male.rows.gt.1000.txt
# gt 100
awk -F ',' '{ for (i=1;i<=NF;i++) if ($i>100) count[i]+=1} END{for (i in count) printf("%d ", count[i]); printf("\n") }' rmats.se.jc.sjc.male.txt > rmats.se.jc.sjc.male.rows.gt.100.txt
awk -F ',' '{ for (i=1;i<=NF;i++) if ($i>100) count[i]+=1} END{for (i in count) printf("%d ", count[i]); printf("\n") }' rmats.se.jc.sjc.female.txt > rmats.se.jc.sjc.female.rows.gt.100.txt
awk -F ',' '{ for (i=1;i<=NF;i++) if ($i>100) count[i]+=1} END{for (i in count) printf("%d ", count[i]); printf("\n") }' rmats.se.jc.ijc.female.txt > rmats.se.jc.ijc.female.rows.gt.100.txt
awk -F ',' '{ for (i=1;i<=NF;i++) if ($i>100) count[i]+=1} END{for (i in count) printf("%d ", count[i]); printf("\n") }' rmats.se.jc.ijc.male.txt > rmats.se.jc.ijc.male.rows.gt.1000.txt
# gt 50
awk -F ',' '{ for (i=1;i<=NF;i++) if ($i>50) count[i]+=1} END{for (i in count) printf("%d ", count[i]); printf("\n") }' rmats.se.jc.sjc.male.txt > rmats.se.jc.sjc.male.rows.gt.50.txt
awk -F ',' '{ for (i=1;i<=NF;i++) if ($i>50) count[i]+=1} END{for (i in count) printf("%d ", count[i]); printf("\n") }' rmats.se.jc.sjc.female.txt > rmats.se.jc.sjc.female.rows.gt.50.txt
awk -F ',' '{ for (i=1;i<=NF;i++) if ($i>50) count[i]+=1} END{for (i in count) printf("%d ", count[i]); printf("\n") }' rmats.se.jc.ijc.female.txt > rmats.se.jc.ijc.female.rows.gt.50.txt
awk -F ',' '{ for (i=1;i<=NF;i++) if ($i>50) count[i]+=1} END{for (i in count) printf("%d ", count[i]); printf("\n") }' rmats.se.jc.ijc.male.txt > rmats.se.jc.ijc.male.rows.gt.50.txt
# gt 10
awk -F ',' '{ for (i=1;i<=NF;i++) if ($i>10) count[i]+=1} END{for (i in count) printf("%d ", count[i]); printf("\n") }' rmats.se.jc.sjc.male.txt > rmats.se.jc.sjc.male.rows.gt.10.txt
awk -F ',' '{ for (i=1;i<=NF;i++) if ($i>10) count[i]+=1} END{for (i in count) printf("%d ", count[i]); printf("\n") }' rmats.se.jc.sjc.female.txt > rmats.se.jc.sjc.female.rows.gt.10.txt
awk -F ',' '{ for (i=1;i<=NF;i++) if ($i>10) count[i]+=1} END{for (i in count) printf("%d ", count[i]); printf("\n") }' rmats.se.jc.ijc.female.txt > rmats.se.jc.ijc.female.rows.gt.10.txt
awk -F ',' '{ for (i=1;i<=NF;i++) if ($i>10) count[i]+=1} END{for (i in count) printf("%d ", count[i]); printf("\n") }' rmats.se.jc.ijc.male.txt > rmats.se.jc.ijc.male.rows.gt.10.txt
# gt 1
awk -F ',' '{ for (i=1;i<=NF;i++) if ($i>1) count[i]+=1} END{for (i in count) printf("%d ", count[i]); printf("\n") }' rmats.se.jc.sjc.male.txt > rmats.se.jc.sjc.male.rows.gt.1.txt
awk -F ',' '{ for (i=1;i<=NF;i++) if ($i>1) count[i]+=1} END{for (i in count) printf("%d ", count[i]); printf("\n") }' rmats.se.jc.sjc.female.txt > rmats.se.jc.sjc.female.rows.gt.1.txt
awk -F ',' '{ for (i=1;i<=NF;i++) if ($i>1) count[i]+=1} END{for (i in count) printf("%d ", count[i]); printf("\n") }' rmats.se.jc.ijc.female.txt > rmats.se.jc.ijc.female.rows.gt.1.txt
awk -F ',' '{ for (i=1;i<=NF;i++) if ($i>1) count[i]+=1} END{for (i in count) printf("%d ", count[i]); printf("\n") }' rmats.se.jc.ijc.male.txt > rmats.se.jc.ijc.male.rows.gt.1.txt

