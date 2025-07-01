# Changelog

All notable changes to this project will be documented here.

## [v1.1] – 2025-06-30
### Added
- Added warning message for positions not present in the reference sequence.
- Added uncertainty scores for reconstructed alleles: UNC values close to 0 present higher confidence of the reconsutrcted ancestral allele, values close to 1 present lower confidence of the ancestral allele, due to mutability and/or lack of data during Ancestral State Reconstruction. 
  
### Fixed
- Corrected column reference to fix script failure when generating output file of multi-sample vcf files. (Thank you, Zhiyong Wang for the bug report!)
- Recognition and annotation of missing data in vcf files

### Changed
- Updated Y-ARS Sequence (Y-ARS v1.1) upon consultation of long-read primate data for ancestral state reconstruction. (There are differences in the alleles compared to version 1.0)

---

## [v1.0] – 2025-04-30
### Added
- Initial release with core functionality.
