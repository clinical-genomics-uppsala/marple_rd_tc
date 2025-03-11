# Changelog

## [0.4.0](https://www.github.com/clinical-genomics-uppsala/marple_rd_tc/compare/v0.3.1...v0.4.0) (2025-03-11)


### Features

* add config with novaseqX exomedepth reference ([34d0941](https://www.github.com/clinical-genomics-uppsala/marple_rd_tc/commit/34d0941efe6ab7b7410894a33e5437b40e7536c8))
* copy the exomedepth config to the results ([a2c8552](https://www.github.com/clinical-genomics-uppsala/marple_rd_tc/commit/a2c85520c3cf0ee8d50c8aae9caaf7c74ac84158))

### [0.3.1](https://www.github.com/clinical-genomics-uppsala/marple_rd_tc/compare/v0.3.0...v0.3.1) (2024-06-04)


### Bug Fixes

* update README.md ([54f54a6](https://www.github.com/clinical-genomics-uppsala/marple_rd_tc/commit/54f54a6a5bdd0c2abd651617cc4e40e609bcc984))

## [0.3.0](https://www.github.com/clinical-genomics-uppsala/marple_rd_tc/compare/v0.2.0...v0.3.0) (2024-04-24)


### Features

* adding PATH_TO_REPO parameter in config ([efad76a](https://www.github.com/clinical-genomics-uppsala/marple_rd_tc/commit/efad76ac1fdd86326bb3d4e76c3e9c3517162808))
* remove SampleSheet.csv for sample order ([3aa4b74](https://www.github.com/clinical-genomics-uppsala/marple_rd_tc/commit/3aa4b748fac13af067e018a8693d946fa18f88cd))


### Bug Fixes

* remove Snakefile in profile ([5dac5ef](https://www.github.com/clinical-genomics-uppsala/marple_rd_tc/commit/5dac5ef8dbc66ac0c70c50fd9e6c689984d375f4))

## [0.2.0](https://www.github.com/clinical-genomics-uppsala/marple_rd_tc/compare/v0.1.4...v0.2.0) (2024-04-04)


### Features

* switch gpu deepvariant to cpu version ([a2a65b0](https://www.github.com/clinical-genomics-uppsala/marple_rd_tc/commit/a2a65b041d0d80d3e2e233fa71a445b79301400e))


### Bug Fixes

* add deepvariant container ([c846250](https://www.github.com/clinical-genomics-uppsala/marple_rd_tc/commit/c84625043cf4b9fe80d36645856b4873b3c2e199))
* add missing tbi files ([665217f](https://www.github.com/clinical-genomics-uppsala/marple_rd_tc/commit/665217f7dbbfd37fa69d12ebe9447d9bfc4c7b13))


### Documentation

* remove gpu from doc add cpu ([3e8e72a](https://www.github.com/clinical-genomics-uppsala/marple_rd_tc/commit/3e8e72a65807a696a18b75ac7b7829de89ab2749))

### [0.1.4](https://www.github.com/clinical-genomics-uppsala/marple_rd_tc/compare/v0.1.3...v0.1.4) (2024-03-14)


### Bug Fixes

* add missing bai to exomedepth, fix typo in start script ([1c58ea7](https://www.github.com/clinical-genomics-uppsala/marple_rd_tc/commit/1c58ea73eb0374aeedd954ad8b1aefaf0d397621))

### [0.1.3](https://www.github.com/clinical-genomics-uppsala/marple_rd_tc/compare/v0.1.2...v0.1.3) (2024-02-13)


### Bug Fixes

* pin pulp to 2.8 or less ([5182cca](https://www.github.com/clinical-genomics-uppsala/marple_rd_tc/commit/5182cca2208cd203b7d35004426f0d76cd06b69f))
* replace total seq (samtools) with total reads (picard) in multiqc ([81eb4a2](https://www.github.com/clinical-genomics-uppsala/marple_rd_tc/commit/81eb4a2fc7eabd93055cbbf098be66a5a274be54))
* run samtools stats without bed, use value in multiqc ([ae637c9](https://www.github.com/clinical-genomics-uppsala/marple_rd_tc/commit/ae637c9ceccc7cc5db699369e4fe5fc6bb94e110))


### Performance Improvements

* use correct num for common docker ([3a2f797](https://www.github.com/clinical-genomics-uppsala/marple_rd_tc/commit/3a2f7978bac89184333481239a684c9225e57daf))


### Documentation

* update samtools stats info ([4ef515f](https://www.github.com/clinical-genomics-uppsala/marple_rd_tc/commit/4ef515f721e6900dcde5dceaa295a28789eb3c51))

### [0.1.2](https://www.github.com/clinical-genomics-uppsala/marple_rd_tc/compare/v0.1.1...v0.1.2) (2023-09-27)


### Bug Fixes

* add config to Results folder ([5b00887](https://www.github.com/clinical-genomics-uppsala/marple_rd_tc/commit/5b00887f7c82691bd483aa35a3a153d5de145ff7))
* update snakemakeprofile and normal panel for exomedepthversion ([f49efe9](https://www.github.com/clinical-genomics-uppsala/marple_rd_tc/commit/f49efe94653afa20d9227bf922c1ac1fc9615353))


### Documentation

* add config to resultsfiles.md ([300fa5b](https://www.github.com/clinical-genomics-uppsala/marple_rd_tc/commit/300fa5b36f5a44d91e98ab1ef28d120592327b5f))
* update mkdocs.yaml ([fa6f819](https://www.github.com/clinical-genomics-uppsala/marple_rd_tc/commit/fa6f819fdb84921b468c18bf57e28af85ee40a21))
* update parabricks and links to main branch of marple ([adf752e](https://www.github.com/clinical-genomics-uppsala/marple_rd_tc/commit/adf752e7bcee51bfd086c4489c7b8bed9cd4220e))
* update README.md ([51bf321](https://www.github.com/clinical-genomics-uppsala/marple_rd_tc/commit/51bf321eb9158360a943fc10b61a4b534e22040e))
* update running.md with some hardware specs ([635e60c](https://www.github.com/clinical-genomics-uppsala/marple_rd_tc/commit/635e60cd4822bdf99bd1cd906c2ce562048fa1ce))
* update some links ([60b4e7a](https://www.github.com/clinical-genomics-uppsala/marple_rd_tc/commit/60b4e7a235c9f6c6a6b547d8bdb887d472a93f8f))

### [0.1.1](https://www.github.com/clinical-genomics-uppsala/marple_rd_tc/compare/v0.1.0...v0.1.1) (2023-09-14)


### Bug Fixes

* remove workflow/rules from mkdocs.yaml ([65c1308](https://www.github.com/clinical-genomics-uppsala/marple_rd_tc/commit/65c13084ec63feaef8c867b2c2bb075298f83120))

## 0.1.0 (2023-09-14)


### Features

* add deepvariant instead of haplotypecaller ([9d8dc06](https://www.github.com/clinical-genomics-uppsala/marple_rd_tc/commit/9d8dc062cc9fd1101200a9bbbbf44444eee9e9c2))
* add g.vcf.gz ([0c2e7e0](https://www.github.com/clinical-genomics-uppsala/marple_rd_tc/commit/0c2e7e004a718631a901d6622aa7e23c38d62daa))
* add pgrs to xlsx file ([9d462d6](https://www.github.com/clinical-genomics-uppsala/marple_rd_tc/commit/9d462d6294e9bd4eef852e483a4fc1b8f234f78b))
* add referece pipeline for exomedepth normal pool ([783fef1](https://www.github.com/clinical-genomics-uppsala/marple_rd_tc/commit/783fef1a846fe2efc1c256cb8c4f818168789fc4))
* first draft pipeline ([f3b3350](https://www.github.com/clinical-genomics-uppsala/marple_rd_tc/commit/f3b3350f4be6e5b8ddbf675042aa4e181b82e329))


### Bug Fixes

* add bedfile to deepvariant ([70fd874](https://www.github.com/clinical-genomics-uppsala/marple_rd_tc/commit/70fd874e384db761c9739166b903049eb1628aeb))
* add exon bedfile to mosdepth and hsmetrics ([d203fbd](https://www.github.com/clinical-genomics-uppsala/marple_rd_tc/commit/d203fbd8d9b1ee333b765055a728fefc5a8d8f1b))
* add missing containers ([31439c5](https://www.github.com/clinical-genomics-uppsala/marple_rd_tc/commit/31439c5e0dc17434c2f2d208c20480052eb471da))
* add reference line to vcf for alissa ([ea01741](https://www.github.com/clinical-genomics-uppsala/marple_rd_tc/commit/ea01741bcf74b0da6475412035fe123bd022a94b))
* make pycodestyle happy ([22de945](https://www.github.com/clinical-genomics-uppsala/marple_rd_tc/commit/22de9452688a5b2a495e74c06e1629b2d4ed9eb7))
* more conda removed ([bc6c0ef](https://www.github.com/clinical-genomics-uppsala/marple_rd_tc/commit/bc6c0ef6867d6da8cc4c5f73c01ccb80a1cf24b0))
* picard hsmetrics on full bed instead of exon only ([4d0b017](https://www.github.com/clinical-genomics-uppsala/marple_rd_tc/commit/4d0b017db4db20b3fde447222b3559613dd88f21))
* remove conda ([d37b9ad](https://www.github.com/clinical-genomics-uppsala/marple_rd_tc/commit/d37b9ad5efbef87a240b459a181aafc60f442b0b))
* remove duplicate rows in low_cov sheet ([179aa49](https://www.github.com/clinical-genomics-uppsala/marple_rd_tc/commit/179aa49349f9bc120e5ed3dd43b6a2f1f3da9a5f))
* remove exomdepth output =1 ([8556e33](https://www.github.com/clinical-genomics-uppsala/marple_rd_tc/commit/8556e33166d968f6e3a07ec50f7af78865c6c80e))
* set coverage as number so table sort works ([2745a78](https://www.github.com/clinical-genomics-uppsala/marple_rd_tc/commit/2745a785d8e71a47c0af5b7747eaa52b69ad1aee))
* spelling misstake pgrs_bed ([9f1bf5e](https://www.github.com/clinical-genomics-uppsala/marple_rd_tc/commit/9f1bf5ee9e7ddc3351349ca132503edeb2df2c2d))
* spelling on requirements file ([8170bdc](https://www.github.com/clinical-genomics-uppsala/marple_rd_tc/commit/8170bdcc22a1279ceb9d815db266688c95a9dd8b))
* update hydra version and update bedtools intersect schema ([00ca936](https://www.github.com/clinical-genomics-uppsala/marple_rd_tc/commit/00ca936d6ec4de55943b177a67d5459a8e0e5e68))
* update parabricks version in config.yaml ([23fc3df](https://www.github.com/clinical-genomics-uppsala/marple_rd_tc/commit/23fc3df0ceb55ceb48711617d00e1a8c3c4c5a85))
* update path to repo ([c1dd22c](https://www.github.com/clinical-genomics-uppsala/marple_rd_tc/commit/c1dd22cdcf2084aaed2e3775d15c7e6fb39cd3dc))


### Documentation

* add changelog to doc hooks.py ([e62441a](https://www.github.com/clinical-genomics-uppsala/marple_rd_tc/commit/e62441aa3840562172e7fcd2098d6bac152896ff))
* add link to github on readthedocs ([db1197e](https://www.github.com/clinical-genomics-uppsala/marple_rd_tc/commit/db1197e75908256a17889b0041e03bcae41eed16))
* add readthedocs ([eeac5f1](https://www.github.com/clinical-genomics-uppsala/marple_rd_tc/commit/eeac5f118f329ace38fbe7afa76583710b0aa474))
* add reference workflow to readthedocs ([739362a](https://www.github.com/clinical-genomics-uppsala/marple_rd_tc/commit/739362a653f5621e9baca5b011be62b8919a5861))
* add rulegraph to hook to cp on build docs ([dadf705](https://www.github.com/clinical-genomics-uppsala/marple_rd_tc/commit/dadf70523847e0c15362cfda23020666aff99742))
* readthedocs add step-by-step ([2913bd3](https://www.github.com/clinical-genomics-uppsala/marple_rd_tc/commit/2913bd325fc76573cd8820aedb85a6c3e08b8285))
* small docs updates ([d2181b8](https://www.github.com/clinical-genomics-uppsala/marple_rd_tc/commit/d2181b8ca58a70d44f74ded9f869846dcca594de))
* structure changes ([704c743](https://www.github.com/clinical-genomics-uppsala/marple_rd_tc/commit/704c743e26d5f04aafb3ff2a4bc544591834973b))
* update CODEOWNERS ([ce8b0a7](https://www.github.com/clinical-genomics-uppsala/marple_rd_tc/commit/ce8b0a70b8445aeb386f3e87afbbaf19cdd58030))
* update hydra readthedocs links ([cc9931f](https://www.github.com/clinical-genomics-uppsala/marple_rd_tc/commit/cc9931fba40b9b39cd0e79f7b7ea8c9b03d83036))
* update paths to git repo ([7fa888e](https://www.github.com/clinical-genomics-uppsala/marple_rd_tc/commit/7fa888eb0974d6a5cfcc8f91f3581cc685fcaa9c))
* update readme ([d8d333c](https://www.github.com/clinical-genomics-uppsala/marple_rd_tc/commit/d8d333c6aff536679a5a9ec6fb5dd9686093b78e))
* update readme with link to readthedocs ([e2a749c](https://www.github.com/clinical-genomics-uppsala/marple_rd_tc/commit/e2a749c5e883da79de1cd17dca9fe3b98f9ea8e2))
* update readthedocs ([ac74ba2](https://www.github.com/clinical-genomics-uppsala/marple_rd_tc/commit/ac74ba2f11e87fd1c70e0141d7e835b7ca46d396))
* update readthedocs to include deepvariant ([ea2f0eb](https://www.github.com/clinical-genomics-uppsala/marple_rd_tc/commit/ea2f0eb180b8c23cf3c056278d55e26d59c102b4))
* update rulegraph ([5ee4577](https://www.github.com/clinical-genomics-uppsala/marple_rd_tc/commit/5ee457745b8f3db43a357259706f3bd6bb383186))
* update running ref ([8f99ae7](https://www.github.com/clinical-genomics-uppsala/marple_rd_tc/commit/8f99ae7f2b802f38617626729ea309de092d3ea5))
* update schemas ([fc25b81](https://www.github.com/clinical-genomics-uppsala/marple_rd_tc/commit/fc25b815590c4b076d2a48b5d6db11b0ca2a5ad5))
