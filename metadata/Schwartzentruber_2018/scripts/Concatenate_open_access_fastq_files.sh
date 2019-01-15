#!/bin/bash
# Concatenate fastq files for open access RNA-seq

mkdir -p concat/SAMEA3234618
cat ERR1775541/ERR1775541_1.fastq.gz ERR1775542/ERR1775542_1.fastq.gz ERR1775543/ERR1775543_1.fastq.gz ERR1775544/ERR1775544_1.fastq.gz > concat/SAMEA3234618/SAMEA3234618_1.fastq.gz
cat ERR1775541/ERR1775541_2.fastq.gz ERR1775542/ERR1775542_2.fastq.gz ERR1775543/ERR1775543_2.fastq.gz ERR1775544/ERR1775544_2.fastq.gz > concat/SAMEA3234618/SAMEA3234618_2.fastq.gz
echo SAMEA3234618

mkdir -p concat/SAMEA3476922
cat ERR1775545/ERR1775545_1.fastq.gz ERR1775551/ERR1775551_1.fastq.gz > concat/SAMEA3476922/SAMEA3476922_1.fastq.gz
cat ERR1775545/ERR1775545_2.fastq.gz ERR1775551/ERR1775551_2.fastq.gz > concat/SAMEA3476922/SAMEA3476922_2.fastq.gz
echo SAMEA3476922

mkdir -p concat/SAMEA3476924
cat ERR1775546/ERR1775546_1.fastq.gz ERR1775552/ERR1775552_1.fastq.gz > concat/SAMEA3476924/SAMEA3476924_1.fastq.gz
cat ERR1775546/ERR1775546_2.fastq.gz ERR1775552/ERR1775552_2.fastq.gz > concat/SAMEA3476924/SAMEA3476924_2.fastq.gz
echo SAMEA3476924

mkdir -p concat/SAMEA3476935
cat ERR1775547/ERR1775547_1.fastq.gz ERR1775553/ERR1775553_1.fastq.gz > concat/SAMEA3476935/SAMEA3476935_1.fastq.gz
cat ERR1775547/ERR1775547_2.fastq.gz ERR1775553/ERR1775553_2.fastq.gz > concat/SAMEA3476935/SAMEA3476935_2.fastq.gz
echo SAMEA3476935

mkdir -p concat/SAMEA3476942
cat ERR1775548/ERR1775548_1.fastq.gz ERR1775554/ERR1775554_1.fastq.gz > concat/SAMEA3476942/SAMEA3476942_1.fastq.gz
cat ERR1775548/ERR1775548_2.fastq.gz ERR1775554/ERR1775554_2.fastq.gz > concat/SAMEA3476942/SAMEA3476942_2.fastq.gz
echo SAMEA3476942

mkdir -p concat/SAMEA3476945
cat ERR1775549/ERR1775549_1.fastq.gz ERR1775555/ERR1775555_1.fastq.gz > concat/SAMEA3476945/SAMEA3476945_1.fastq.gz
cat ERR1775549/ERR1775549_2.fastq.gz ERR1775555/ERR1775555_2.fastq.gz > concat/SAMEA3476945/SAMEA3476945_2.fastq.gz
echo SAMEA3476945

mkdir -p concat/SAMEA3476950
cat ERR1775550/ERR1775550_1.fastq.gz ERR1775556/ERR1775556_1.fastq.gz > concat/SAMEA3476950/SAMEA3476950_1.fastq.gz
cat ERR1775550/ERR1775550_2.fastq.gz ERR1775556/ERR1775556_2.fastq.gz > concat/SAMEA3476950/SAMEA3476950_2.fastq.gz
echo SAMEA3476950

mkdir -p concat/SAMEA3476967
cat ERR1775557/ERR1775557_1.fastq.gz ERR1775562/ERR1775562_1.fastq.gz > concat/SAMEA3476967/SAMEA3476967_1.fastq.gz
cat ERR1775557/ERR1775557_2.fastq.gz ERR1775562/ERR1775562_2.fastq.gz > concat/SAMEA3476967/SAMEA3476967_2.fastq.gz
echo SAMEA3476967

mkdir -p concat/SAMEA3476969
cat ERR1775558/ERR1775558_1.fastq.gz ERR1775563/ERR1775563_1.fastq.gz > concat/SAMEA3476969/SAMEA3476969_1.fastq.gz
cat ERR1775558/ERR1775558_2.fastq.gz ERR1775563/ERR1775563_2.fastq.gz > concat/SAMEA3476969/SAMEA3476969_2.fastq.gz
echo SAMEA3476969

mkdir -p concat/SAMEA3476974
cat ERR1775559/ERR1775559_1.fastq.gz ERR1775564/ERR1775564_1.fastq.gz > concat/SAMEA3476974/SAMEA3476974_1.fastq.gz
cat ERR1775559/ERR1775559_2.fastq.gz ERR1775564/ERR1775564_2.fastq.gz > concat/SAMEA3476974/SAMEA3476974_2.fastq.gz
echo SAMEA3476974

mkdir -p concat/SAMEA3476976
cat ERR1775560/ERR1775560_1.fastq.gz ERR1775565/ERR1775565_1.fastq.gz > concat/SAMEA3476976/SAMEA3476976_1.fastq.gz
cat ERR1775560/ERR1775560_2.fastq.gz ERR1775565/ERR1775565_2.fastq.gz > concat/SAMEA3476976/SAMEA3476976_2.fastq.gz
echo SAMEA3476976

mkdir -p concat/SAMEA3476984
cat ERR1775561/ERR1775561_1.fastq.gz ERR1775566/ERR1775566_1.fastq.gz > concat/SAMEA3476984/SAMEA3476984_1.fastq.gz
cat ERR1775561/ERR1775561_2.fastq.gz ERR1775566/ERR1775566_2.fastq.gz > concat/SAMEA3476984/SAMEA3476984_2.fastq.gz
echo SAMEA3476984

mkdir -p concat/SAMEA3864975
cat ERR1775567/ERR1775567_1.fastq.gz ERR1775579/ERR1775579_1.fastq.gz ERR1775591/ERR1775591_1.fastq.gz > concat/SAMEA3864975/SAMEA3864975_1.fastq.gz
cat ERR1775567/ERR1775567_2.fastq.gz ERR1775579/ERR1775579_2.fastq.gz ERR1775591/ERR1775591_2.fastq.gz > concat/SAMEA3864975/SAMEA3864975_2.fastq.gz
echo SAMEA3864975

mkdir -p concat/SAMEA3864976
cat ERR1775568/ERR1775568_1.fastq.gz ERR1775580/ERR1775580_1.fastq.gz ERR1775592/ERR1775592_1.fastq.gz > concat/SAMEA3864976/SAMEA3864976_1.fastq.gz
cat ERR1775568/ERR1775568_2.fastq.gz ERR1775580/ERR1775580_2.fastq.gz ERR1775592/ERR1775592_2.fastq.gz > concat/SAMEA3864976/SAMEA3864976_2.fastq.gz
echo SAMEA3864976

mkdir -p concat/SAMEA3864978
cat ERR1775569/ERR1775569_1.fastq.gz ERR1775581/ERR1775581_1.fastq.gz ERR1775593/ERR1775593_1.fastq.gz > concat/SAMEA3864978/SAMEA3864978_1.fastq.gz
cat ERR1775569/ERR1775569_2.fastq.gz ERR1775581/ERR1775581_2.fastq.gz ERR1775593/ERR1775593_2.fastq.gz > concat/SAMEA3864978/SAMEA3864978_2.fastq.gz
echo SAMEA3864978

mkdir -p concat/SAMEA3864979
cat ERR1775570/ERR1775570_1.fastq.gz ERR1775582/ERR1775582_1.fastq.gz ERR1775594/ERR1775594_1.fastq.gz > concat/SAMEA3864979/SAMEA3864979_1.fastq.gz
cat ERR1775570/ERR1775570_2.fastq.gz ERR1775582/ERR1775582_2.fastq.gz ERR1775594/ERR1775594_2.fastq.gz > concat/SAMEA3864979/SAMEA3864979_2.fastq.gz
echo SAMEA3864979

mkdir -p concat/SAMEA3864980
cat ERR1775571/ERR1775571_1.fastq.gz ERR1775583/ERR1775583_1.fastq.gz ERR1775595/ERR1775595_1.fastq.gz > concat/SAMEA3864980/SAMEA3864980_1.fastq.gz
cat ERR1775571/ERR1775571_2.fastq.gz ERR1775583/ERR1775583_2.fastq.gz ERR1775595/ERR1775595_2.fastq.gz > concat/SAMEA3864980/SAMEA3864980_2.fastq.gz
echo SAMEA3864980

mkdir -p concat/SAMEA3864981
cat ERR1775572/ERR1775572_1.fastq.gz ERR1775584/ERR1775584_1.fastq.gz ERR1775596/ERR1775596_1.fastq.gz > concat/SAMEA3864981/SAMEA3864981_1.fastq.gz
cat ERR1775572/ERR1775572_2.fastq.gz ERR1775584/ERR1775584_2.fastq.gz ERR1775596/ERR1775596_2.fastq.gz > concat/SAMEA3864981/SAMEA3864981_2.fastq.gz
echo SAMEA3864981

mkdir -p concat/SAMEA3864982
cat ERR1775573/ERR1775573_1.fastq.gz ERR1775585/ERR1775585_1.fastq.gz ERR1775597/ERR1775597_1.fastq.gz > concat/SAMEA3864982/SAMEA3864982_1.fastq.gz
cat ERR1775573/ERR1775573_2.fastq.gz ERR1775585/ERR1775585_2.fastq.gz ERR1775597/ERR1775597_2.fastq.gz > concat/SAMEA3864982/SAMEA3864982_2.fastq.gz
echo SAMEA3864982

mkdir -p concat/SAMEA3864969
cat ERR1775574/ERR1775574_1.fastq.gz ERR1775586/ERR1775586_1.fastq.gz ERR1775598/ERR1775598_1.fastq.gz > concat/SAMEA3864969/SAMEA3864969_1.fastq.gz
cat ERR1775574/ERR1775574_2.fastq.gz ERR1775586/ERR1775586_2.fastq.gz ERR1775598/ERR1775598_2.fastq.gz > concat/SAMEA3864969/SAMEA3864969_2.fastq.gz
echo SAMEA3864969

mkdir -p concat/SAMEA3864970
cat ERR1775575/ERR1775575_1.fastq.gz ERR1775587/ERR1775587_1.fastq.gz ERR1775599/ERR1775599_1.fastq.gz > concat/SAMEA3864970/SAMEA3864970_1.fastq.gz
cat ERR1775575/ERR1775575_2.fastq.gz ERR1775587/ERR1775587_2.fastq.gz ERR1775599/ERR1775599_2.fastq.gz > concat/SAMEA3864970/SAMEA3864970_2.fastq.gz
echo SAMEA3864970

mkdir -p concat/SAMEA3864971
cat ERR1775576/ERR1775576_1.fastq.gz ERR1775588/ERR1775588_1.fastq.gz ERR1775600/ERR1775600_1.fastq.gz > concat/SAMEA3864971/SAMEA3864971_1.fastq.gz
cat ERR1775576/ERR1775576_2.fastq.gz ERR1775588/ERR1775588_2.fastq.gz ERR1775600/ERR1775600_2.fastq.gz > concat/SAMEA3864971/SAMEA3864971_2.fastq.gz
echo SAMEA3864971

mkdir -p concat/SAMEA3864972
cat ERR1775577/ERR1775577_1.fastq.gz ERR1775589/ERR1775589_1.fastq.gz ERR1775601/ERR1775601_1.fastq.gz > concat/SAMEA3864972/SAMEA3864972_1.fastq.gz
cat ERR1775577/ERR1775577_2.fastq.gz ERR1775589/ERR1775589_2.fastq.gz ERR1775601/ERR1775601_2.fastq.gz > concat/SAMEA3864972/SAMEA3864972_2.fastq.gz
echo SAMEA3864972

mkdir -p concat/SAMEA3864973
cat ERR1775578/ERR1775578_1.fastq.gz ERR1775590/ERR1775590_1.fastq.gz ERR1775602/ERR1775602_1.fastq.gz > concat/SAMEA3864973/SAMEA3864973_1.fastq.gz
cat ERR1775578/ERR1775578_2.fastq.gz ERR1775590/ERR1775590_2.fastq.gz ERR1775602/ERR1775602_2.fastq.gz > concat/SAMEA3864973/SAMEA3864973_2.fastq.gz
echo SAMEA3864973

mkdir -p concat/SAMEA3864991
cat ERR1775603/ERR1775603_1.fastq.gz ERR1775617/ERR1775617_1.fastq.gz ERR1775631/ERR1775631_1.fastq.gz > concat/SAMEA3864991/SAMEA3864991_1.fastq.gz
cat ERR1775603/ERR1775603_2.fastq.gz ERR1775617/ERR1775617_2.fastq.gz ERR1775631/ERR1775631_2.fastq.gz > concat/SAMEA3864991/SAMEA3864991_2.fastq.gz
echo SAMEA3864991

mkdir -p concat/SAMEA3864992
cat ERR1775604/ERR1775604_1.fastq.gz ERR1775618/ERR1775618_1.fastq.gz ERR1775632/ERR1775632_1.fastq.gz > concat/SAMEA3864992/SAMEA3864992_1.fastq.gz
cat ERR1775604/ERR1775604_2.fastq.gz ERR1775618/ERR1775618_2.fastq.gz ERR1775632/ERR1775632_2.fastq.gz > concat/SAMEA3864992/SAMEA3864992_2.fastq.gz
echo SAMEA3864992

mkdir -p concat/SAMEA3864993
cat ERR1775605/ERR1775605_1.fastq.gz ERR1775619/ERR1775619_1.fastq.gz ERR1775633/ERR1775633_1.fastq.gz > concat/SAMEA3864993/SAMEA3864993_1.fastq.gz
cat ERR1775605/ERR1775605_2.fastq.gz ERR1775619/ERR1775619_2.fastq.gz ERR1775633/ERR1775633_2.fastq.gz > concat/SAMEA3864993/SAMEA3864993_2.fastq.gz
echo SAMEA3864993

mkdir -p concat/SAMEA3864994
cat ERR1775606/ERR1775606_1.fastq.gz ERR1775620/ERR1775620_1.fastq.gz ERR1775634/ERR1775634_1.fastq.gz > concat/SAMEA3864994/SAMEA3864994_1.fastq.gz
cat ERR1775606/ERR1775606_2.fastq.gz ERR1775620/ERR1775620_2.fastq.gz ERR1775634/ERR1775634_2.fastq.gz > concat/SAMEA3864994/SAMEA3864994_2.fastq.gz
echo SAMEA3864994

mkdir -p concat/SAMEA3864995
cat ERR1775607/ERR1775607_1.fastq.gz ERR1775621/ERR1775621_1.fastq.gz ERR1775635/ERR1775635_1.fastq.gz > concat/SAMEA3864995/SAMEA3864995_1.fastq.gz
cat ERR1775607/ERR1775607_2.fastq.gz ERR1775621/ERR1775621_2.fastq.gz ERR1775635/ERR1775635_2.fastq.gz > concat/SAMEA3864995/SAMEA3864995_2.fastq.gz
echo SAMEA3864995

mkdir -p concat/SAMEA3864996
cat ERR1775608/ERR1775608_1.fastq.gz ERR1775622/ERR1775622_1.fastq.gz ERR1775636/ERR1775636_1.fastq.gz > concat/SAMEA3864996/SAMEA3864996_1.fastq.gz
cat ERR1775608/ERR1775608_2.fastq.gz ERR1775622/ERR1775622_2.fastq.gz ERR1775636/ERR1775636_2.fastq.gz > concat/SAMEA3864996/SAMEA3864996_2.fastq.gz
echo SAMEA3864996

mkdir -p concat/SAMEA3864998
cat ERR1775609/ERR1775609_1.fastq.gz ERR1775623/ERR1775623_1.fastq.gz ERR1775637/ERR1775637_1.fastq.gz > concat/SAMEA3864998/SAMEA3864998_1.fastq.gz
cat ERR1775609/ERR1775609_2.fastq.gz ERR1775623/ERR1775623_2.fastq.gz ERR1775637/ERR1775637_2.fastq.gz > concat/SAMEA3864998/SAMEA3864998_2.fastq.gz
echo SAMEA3864998

mkdir -p concat/SAMEA3864984
cat ERR1775610/ERR1775610_1.fastq.gz ERR1775624/ERR1775624_1.fastq.gz ERR1775638/ERR1775638_1.fastq.gz > concat/SAMEA3864984/SAMEA3864984_1.fastq.gz
cat ERR1775610/ERR1775610_2.fastq.gz ERR1775624/ERR1775624_2.fastq.gz ERR1775638/ERR1775638_2.fastq.gz > concat/SAMEA3864984/SAMEA3864984_2.fastq.gz
echo SAMEA3864984

mkdir -p concat/SAMEA3864985
cat ERR1775611/ERR1775611_1.fastq.gz ERR1775625/ERR1775625_1.fastq.gz ERR1775639/ERR1775639_1.fastq.gz > concat/SAMEA3864985/SAMEA3864985_1.fastq.gz
cat ERR1775611/ERR1775611_2.fastq.gz ERR1775625/ERR1775625_2.fastq.gz ERR1775639/ERR1775639_2.fastq.gz > concat/SAMEA3864985/SAMEA3864985_2.fastq.gz
echo SAMEA3864985

mkdir -p concat/SAMEA3864986
cat ERR1775612/ERR1775612_1.fastq.gz ERR1775626/ERR1775626_1.fastq.gz ERR1775640/ERR1775640_1.fastq.gz > concat/SAMEA3864986/SAMEA3864986_1.fastq.gz
cat ERR1775612/ERR1775612_2.fastq.gz ERR1775626/ERR1775626_2.fastq.gz ERR1775640/ERR1775640_2.fastq.gz > concat/SAMEA3864986/SAMEA3864986_2.fastq.gz
echo SAMEA3864986

mkdir -p concat/SAMEA3864987
cat ERR1775613/ERR1775613_1.fastq.gz ERR1775627/ERR1775627_1.fastq.gz ERR1775641/ERR1775641_1.fastq.gz > concat/SAMEA3864987/SAMEA3864987_1.fastq.gz
cat ERR1775613/ERR1775613_2.fastq.gz ERR1775627/ERR1775627_2.fastq.gz ERR1775641/ERR1775641_2.fastq.gz > concat/SAMEA3864987/SAMEA3864987_2.fastq.gz
echo SAMEA3864987

mkdir -p concat/SAMEA3864988
cat ERR1775614/ERR1775614_1.fastq.gz ERR1775628/ERR1775628_1.fastq.gz ERR1775642/ERR1775642_1.fastq.gz > concat/SAMEA3864988/SAMEA3864988_1.fastq.gz
cat ERR1775614/ERR1775614_2.fastq.gz ERR1775628/ERR1775628_2.fastq.gz ERR1775642/ERR1775642_2.fastq.gz > concat/SAMEA3864988/SAMEA3864988_2.fastq.gz
echo SAMEA3864988

mkdir -p concat/SAMEA3864989
cat ERR1775615/ERR1775615_1.fastq.gz ERR1775629/ERR1775629_1.fastq.gz ERR1775643/ERR1775643_1.fastq.gz > concat/SAMEA3864989/SAMEA3864989_1.fastq.gz
cat ERR1775615/ERR1775615_2.fastq.gz ERR1775629/ERR1775629_2.fastq.gz ERR1775643/ERR1775643_2.fastq.gz > concat/SAMEA3864989/SAMEA3864989_2.fastq.gz
echo SAMEA3864989

mkdir -p concat/SAMEA3864990
cat ERR1775616/ERR1775616_1.fastq.gz ERR1775630/ERR1775630_1.fastq.gz ERR1775644/ERR1775644_1.fastq.gz > concat/SAMEA3864990/SAMEA3864990_1.fastq.gz
cat ERR1775616/ERR1775616_2.fastq.gz ERR1775630/ERR1775630_2.fastq.gz ERR1775644/ERR1775644_2.fastq.gz > concat/SAMEA3864990/SAMEA3864990_2.fastq.gz
echo SAMEA3864990

mkdir -p concat/SAMEA3717998
cat ERR1775645/ERR1775645_1.fastq.gz ERR1775658/ERR1775658_1.fastq.gz ERR1775671/ERR1775671_1.fastq.gz ERR1775684/ERR1775684_1.fastq.gz > concat/SAMEA3717998/SAMEA3717998_1.fastq.gz
cat ERR1775645/ERR1775645_2.fastq.gz ERR1775658/ERR1775658_2.fastq.gz ERR1775671/ERR1775671_2.fastq.gz ERR1775684/ERR1775684_2.fastq.gz > concat/SAMEA3717998/SAMEA3717998_2.fastq.gz
echo SAMEA3717998

mkdir -p concat/SAMEA3718000
cat ERR1775646/ERR1775646_1.fastq.gz ERR1775659/ERR1775659_1.fastq.gz ERR1775672/ERR1775672_1.fastq.gz ERR1775685/ERR1775685_1.fastq.gz > concat/SAMEA3718000/SAMEA3718000_1.fastq.gz
cat ERR1775646/ERR1775646_2.fastq.gz ERR1775659/ERR1775659_2.fastq.gz ERR1775672/ERR1775672_2.fastq.gz ERR1775685/ERR1775685_2.fastq.gz > concat/SAMEA3718000/SAMEA3718000_2.fastq.gz
echo SAMEA3718000

mkdir -p concat/SAMEA3718001
cat ERR1775647/ERR1775647_1.fastq.gz ERR1775660/ERR1775660_1.fastq.gz ERR1775673/ERR1775673_1.fastq.gz ERR1775686/ERR1775686_1.fastq.gz > concat/SAMEA3718001/SAMEA3718001_1.fastq.gz
cat ERR1775647/ERR1775647_2.fastq.gz ERR1775660/ERR1775660_2.fastq.gz ERR1775673/ERR1775673_2.fastq.gz ERR1775686/ERR1775686_2.fastq.gz > concat/SAMEA3718001/SAMEA3718001_2.fastq.gz
echo SAMEA3718001

mkdir -p concat/SAMEA3718002
cat ERR1775648/ERR1775648_1.fastq.gz ERR1775661/ERR1775661_1.fastq.gz ERR1775674/ERR1775674_1.fastq.gz ERR1775687/ERR1775687_1.fastq.gz > concat/SAMEA3718002/SAMEA3718002_1.fastq.gz
cat ERR1775648/ERR1775648_2.fastq.gz ERR1775661/ERR1775661_2.fastq.gz ERR1775674/ERR1775674_2.fastq.gz ERR1775687/ERR1775687_2.fastq.gz > concat/SAMEA3718002/SAMEA3718002_2.fastq.gz
echo SAMEA3718002

mkdir -p concat/SAMEA3718653
cat ERR1775649/ERR1775649_1.fastq.gz ERR1775662/ERR1775662_1.fastq.gz ERR1775675/ERR1775675_1.fastq.gz ERR1775688/ERR1775688_1.fastq.gz > concat/SAMEA3718653/SAMEA3718653_1.fastq.gz
cat ERR1775649/ERR1775649_2.fastq.gz ERR1775662/ERR1775662_2.fastq.gz ERR1775675/ERR1775675_2.fastq.gz ERR1775688/ERR1775688_2.fastq.gz > concat/SAMEA3718653/SAMEA3718653_2.fastq.gz
echo SAMEA3718653

mkdir -p concat/SAMEA3718658
cat ERR1775650/ERR1775650_1.fastq.gz ERR1775663/ERR1775663_1.fastq.gz ERR1775676/ERR1775676_1.fastq.gz ERR1775689/ERR1775689_1.fastq.gz > concat/SAMEA3718658/SAMEA3718658_1.fastq.gz
cat ERR1775650/ERR1775650_2.fastq.gz ERR1775663/ERR1775663_2.fastq.gz ERR1775676/ERR1775676_2.fastq.gz ERR1775689/ERR1775689_2.fastq.gz > concat/SAMEA3718658/SAMEA3718658_2.fastq.gz
echo SAMEA3718658

mkdir -p concat/SAMEA3718654
cat ERR1775651/ERR1775651_1.fastq.gz ERR1775664/ERR1775664_1.fastq.gz ERR1775677/ERR1775677_1.fastq.gz ERR1775690/ERR1775690_1.fastq.gz > concat/SAMEA3718654/SAMEA3718654_1.fastq.gz
cat ERR1775651/ERR1775651_2.fastq.gz ERR1775664/ERR1775664_2.fastq.gz ERR1775677/ERR1775677_2.fastq.gz ERR1775690/ERR1775690_2.fastq.gz > concat/SAMEA3718654/SAMEA3718654_2.fastq.gz
echo SAMEA3718654

mkdir -p concat/SAMEA3718656
cat ERR1775652/ERR1775652_1.fastq.gz ERR1775665/ERR1775665_1.fastq.gz ERR1775678/ERR1775678_1.fastq.gz ERR1775691/ERR1775691_1.fastq.gz > concat/SAMEA3718656/SAMEA3718656_1.fastq.gz
cat ERR1775652/ERR1775652_2.fastq.gz ERR1775665/ERR1775665_2.fastq.gz ERR1775678/ERR1775678_2.fastq.gz ERR1775691/ERR1775691_2.fastq.gz > concat/SAMEA3718656/SAMEA3718656_2.fastq.gz
echo SAMEA3718656

mkdir -p concat/SAMEA3718661
cat ERR1775653/ERR1775653_1.fastq.gz ERR1775666/ERR1775666_1.fastq.gz ERR1775679/ERR1775679_1.fastq.gz ERR1775692/ERR1775692_1.fastq.gz > concat/SAMEA3718661/SAMEA3718661_1.fastq.gz
cat ERR1775653/ERR1775653_2.fastq.gz ERR1775666/ERR1775666_2.fastq.gz ERR1775679/ERR1775679_2.fastq.gz ERR1775692/ERR1775692_2.fastq.gz > concat/SAMEA3718661/SAMEA3718661_2.fastq.gz
echo SAMEA3718661

mkdir -p concat/SAMEA3718662
cat ERR1775654/ERR1775654_1.fastq.gz ERR1775667/ERR1775667_1.fastq.gz ERR1775680/ERR1775680_1.fastq.gz ERR1775693/ERR1775693_1.fastq.gz > concat/SAMEA3718662/SAMEA3718662_1.fastq.gz
cat ERR1775654/ERR1775654_2.fastq.gz ERR1775667/ERR1775667_2.fastq.gz ERR1775680/ERR1775680_2.fastq.gz ERR1775693/ERR1775693_2.fastq.gz > concat/SAMEA3718662/SAMEA3718662_2.fastq.gz
echo SAMEA3718662

mkdir -p concat/SAMEA3718665
cat ERR1775655/ERR1775655_1.fastq.gz ERR1775668/ERR1775668_1.fastq.gz ERR1775681/ERR1775681_1.fastq.gz ERR1775694/ERR1775694_1.fastq.gz > concat/SAMEA3718665/SAMEA3718665_1.fastq.gz
cat ERR1775655/ERR1775655_2.fastq.gz ERR1775668/ERR1775668_2.fastq.gz ERR1775681/ERR1775681_2.fastq.gz ERR1775694/ERR1775694_2.fastq.gz > concat/SAMEA3718665/SAMEA3718665_2.fastq.gz
echo SAMEA3718665

mkdir -p concat/SAMEA3718666
cat ERR1775656/ERR1775656_1.fastq.gz ERR1775669/ERR1775669_1.fastq.gz ERR1775682/ERR1775682_1.fastq.gz ERR1775695/ERR1775695_1.fastq.gz > concat/SAMEA3718666/SAMEA3718666_1.fastq.gz
cat ERR1775656/ERR1775656_2.fastq.gz ERR1775669/ERR1775669_2.fastq.gz ERR1775682/ERR1775682_2.fastq.gz ERR1775695/ERR1775695_2.fastq.gz > concat/SAMEA3718666/SAMEA3718666_2.fastq.gz
echo SAMEA3718666

mkdir -p concat/SAMEA3718672
cat ERR1775657/ERR1775657_1.fastq.gz ERR1775670/ERR1775670_1.fastq.gz ERR1775683/ERR1775683_1.fastq.gz ERR1775696/ERR1775696_1.fastq.gz > concat/SAMEA3718672/SAMEA3718672_1.fastq.gz
cat ERR1775657/ERR1775657_2.fastq.gz ERR1775670/ERR1775670_2.fastq.gz ERR1775683/ERR1775683_2.fastq.gz ERR1775696/ERR1775696_2.fastq.gz > concat/SAMEA3718672/SAMEA3718672_2.fastq.gz
echo SAMEA3718672
