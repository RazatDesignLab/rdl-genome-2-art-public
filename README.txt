Updated December 2023

This is a repository within the ACTIV TRACE initiative that houses a comprehensive collection of datasets related to SARS-CoV-2.
The processing of SARS-CoV-2 Sequence Read Archive (SRA) files has been optimized to identify genetic variations in viral samples.
This information is then presented in the Variant Call Format (VCF). Each VCF file corresponds to the SRA parent-run accessionID.
Additionally, the metadata is available in the parquet format, making it easier to search and filter using the AWS Athena Service.
The SARS-CoV-2 Variant Calling Pipeline is designed to handle new data every six hours, with updates to the AWS ODP bucket occurring daily.

The "vcf/" directory contains VCF files generated from SRA data and organized based on the parent-run accessionID.
For detailed information on how VCF files are generated, you can refer to this link: https://www.ncbi.nlm.nih.gov/sra/docs/sars-cov-2-variant-calling/.

To access metadata for this dataset, you can find it using the command: "aws s3 ls s3://sra-pub-sars-cov2-metadata-us-east-1 --no-sign-request".

The "sra-src/" and "run/" directories contained the normalized data and original format datasets.
Starting April 2023, normalized data and source files can be accessed from the NIH NCBI SRA on AWS registry (https://registry.opendata.aws/ncbi-sra/).
Original format data will be moved to AWS cold storage.
You can retrieve original format data using the Cloud Data Delivery Service (https://www.ncbi.nlm.nih.gov/sra/docs/data-delivery/).