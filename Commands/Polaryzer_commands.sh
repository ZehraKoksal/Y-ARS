#In bash shell, Polaryzer v1.0 was run using the following commands:

#For vcf files aligned to GRCh37:
python polarYzer.py -chromosome NC_000024.9 -reference GRCh37 -input_single_vcf ./37/

#For vcf files aligned to GRCh38:
python polarYzer.py -chromosome NC_000024.10 -reference GRCh37 -input_single_vcf ./38/

#For vcf files aligned to T2T:
python polarYzer.py -chromosome NC_060948.1 -reference T2T -input_single_vcf ./T2T/

#For vcf files aligned to Y-ARS: #Reference T2T can be selected, because Y-ARS uses the same coordinate system!
python polarYzer.py -chromosome NC_060948.1 -reference T2T -input_single_vcf ./YARS/




#For vcf files aligned to GRCh37 and THEN lifted over to T2T:
python polarYzer.py -chromosome NC_000024.9 -reference T2T -input_single_vcf ./t2t_lifted/37/

#For vcf files aligned to GRCh38 and THEN lifted over to T2T:
python polarYzer.py -chromosome NC_000024.10 -reference T2T -input_single_vcf ./t2t_lifted/38/
