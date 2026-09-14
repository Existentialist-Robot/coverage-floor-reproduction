# MANIFEST

Every file in this repository, with one line of purpose and its SHA-256.

## This bundle supersedes v1.0.0 (DOI 10.5281/zenodo.22699236)

Version 1.0.0 of this repository, archived as
[10.5281/zenodo.22699236](https://doi.org/10.5281/zenodo.22699236), **predates the corrected**
**analyses.** After a pre-submission methods review, the manuscript was revised on 2026-09-11 and
the inferential record was recomputed: the post-hoc comparison family is now enumerated (64
comparisons) and corrected for multiplicity under Holm and Benjamini-Hochberg, and several rates
and intervals that had been computed over correlated runs were rebuilt over independent seed
groups. Conclusions that did not survive were withdrawn rather than softened.

A reviewer who follows the DOI lands on superseded material. The Zenodo record is deliberately
left unchanged rather than overwritten; a new archived version will be deposited separately.
**Until that deposit exists, this repository is the current reproduction material.**

The pre-review record is retained verbatim at
`sim/history/RESULTS-before-review-2026-09-11.md`, so the correction can be read rather than
taken on trust.

## Verification

```bash
python sim/review_reanalysis.py --check
python sim/recompute_group_fpr.py
python sim/emit_missingness_arm_interval.py --check
python sim/test_review_reanalysis.py
python sim/verify_reported_numbers.py
```

All five exit zero on the files hashed below. The last is the legacy verifier: 3,873 checks
binding 128 historical ledger claims. It binds the legacy record, not the revised manuscript.

## Scope of this manifest

The table lists every file in the repository except `MANIFEST.md` and `MANIFEST.json`
themselves. `MANIFEST.json` carries the same hashes in machine-readable form.
All data are synthetic simulation output; the bundle contains no participant records.

| File | Purpose | Bytes | SHA-256 |
|---|---|---:|---|
| `.gitignore` | Git ignore rules for this bundle. | 19 | `862263fa1f46c20f0d1e4dac5ffcc75abd55c08211b2c3864c5f8764b9d87793` |
| `.zenodo.json` | Zenodo deposit metadata for the archived release (authors, title, licences). Unchanged from v1.0.0; no new DOI is minted by this repository. | 696 | `abbd0186035343890c4fcef0f4bdf8651bee7ab3aac9fdf5c2536cd0330066e9` |
| `LICENSE` | MIT licence, covering the Python code in this bundle. | 1103 | `c01d0ff537678069d7abd43836cc3f89dcfe31c25e38d8d4937fc0b139cf71d3` |
| `LICENSE-DATA` | Creative Commons Attribution 4.0 licence, covering the aggregates, configs, results record and figures. | 1134 | `466b3149874734fd9d8bfa830fcedb15d079180e471a2bed0c7cb68a42a305da` |
| `README.md` | Bundle overview: what supersedes v1.0.0, what changed, and how to run the revised reanalysis. | 9318 | `8c8ad29a02ffb5e4ebc838ea8c143943f4fca1164ae53766bf36ba8c5110a02e` |
| `sim/README.md` | Campaign documentation: what each script is, and the order a full re-run executes in. | 6524 | `fbffc42c5ba9681d6166de23ebe7ab24f6091a19db53a56afc3839fc61319afd` |
| `sim/REPRODUCE.md` | The reproduction route for the revised analysis: requirements, commands and statistical scope. | 1282 | `ca53a65e7b1eafca294d11ae15189a494a376d168ebf07aa7ebdaffbbc47f136` |
| `sim/RESULTS.md` | The current results record: revised claims, sample units, exclusions, caveats and negative results. | 14460 | `16efe9fffbb000b4ff9e93c4a34900787da70aec6e7737a436d92f9174d51e64` |
| `sim/aggregates/aux_arms.json` | Committed run aggregate: auxiliary revocation and propensity-correction arms, campaign-1 rule. | 242927 | `eb32cc973839fef929ac1869b609bc7bd92bb4ba0a1abd18cc2910c38b6f2d83` |
| `sim/aggregates/aux_arms_recal.json` | Committed run aggregate: auxiliary arms rescored under the recalibrated rule. | 244007 | `0e094af6f717d68dc90282edaff69035bb138668bb63755295501f79ee247b5d` |
| `sim/aggregates/calibration_null_bothclass.json` | Committed run aggregate: dedicated null scored under the both-class rule variant. | 103388 | `96533442d3a343b7409fec5f3ec9e6865e157af02f6ea865c56faae2c86ec937` |
| `sim/aggregates/calibration_null_large.json` | Committed run aggregate: the 600-run dedicated null used for recalibration. | 232441 | `268068044fad1850cd4906dd4a1cbd0037202801b3d53c4b1aab2b0e0285dcdf` |
| `sim/aggregates/calibration_probe.json` | Committed run aggregate: exploratory calibration probe. | 2937 | `9b3c00d9bed4b97dd3ee7c204345ef3ea0e16ba5a416866f1830de001c7bf491` |
| `sim/aggregates/calibration_rule_grids.json` | Committed run aggregate: threshold grids searched during calibration. | 18993 | `8961873c06345f0fbe170a3477d41e5e45ca8e59622373860cc8c66012a8fb25` |
| `sim/aggregates/cluster_lead_20260906.json` | Derived aggregate: cluster-level lead-time analysis (2026-09-06). | 4232 | `11b72072a3450dd28db90dbede31042c863c4f8d6d6b02513a63d9727599a63e` |
| `sim/aggregates/cluster_matching_20260904.json` | Derived aggregate: group construction and matching for the cluster analyses (2026-09-04). | 16411 | `fe542b45b7105336d7af363cee0a948abb39627c0e6314e057ecc520a748724e` |
| `sim/aggregates/cluster_recheck_20260904.json` | Derived aggregate: independent recheck of the cluster construction (2026-09-04). | 3518 | `c79a9563d1181c9877d3c27607d8e254a0e5e9e142d3bdc29dc2b9f30f918b68` |
| `sim/aggregates/coverage_finegrid_block_b_20260904.json` | Committed run aggregate: fine coverage grid, block B (2026-09-04). | 12086963 | `5a506fceb66373d298ef1abac70eef8b526562d3162840a338719031f70faee9` |
| `sim/aggregates/ea_headline.json` | Committed run aggregate: headline coverage-and-noise campaign, campaign-1 rule. | 1786154 | `9843cf66074965335efd6ae9cb74e01a78ebf38c793c1282684e9114ebffa408` |
| `sim/aggregates/ea_headline_recal.json` | Committed run aggregate: headline campaign rescored under the recalibrated rule. The source for every headline number the paper reports. | 2986379 | `0c9b1b53e4ae6d5de5d285586116c7d5aef5e0bb15e402204fa2b4c29f169233` |
| `sim/aggregates/eb_joint_sweep.json` | Committed run aggregate: joint coverage-by-coupling sweep. | 894509 | `fe920de35cc3698198c98ffce67a61c735093f878a117dd0e25a73ce40e4d252` |
| `sim/aggregates/eb_joint_sweep_n500.json` | Committed run aggregate: joint coverage-by-coupling sweep at 500 matched pairs per cell. | 3695293 | `3cab4105094ad91efd10acf91a5e6e208b8cff948700ef6a339e0bd334a554be` |
| `sim/aggregates/gate1_bifurcation.json` | Committed run aggregate: bifurcation-validation results (saddle-node, bimodality, hysteresis, relaxation time). | 239225 | `5f3e4e71c706ab31e2715fc0d84f3ebafb195b5542412c070f51a724261962d3` |
| `sim/aggregates/graph_detector_calibration.json` | Committed run aggregate: graph-detector calibration on a dedicated null. | 228018 | `82d42fa35b76da4a0fab352cdf6c97e1228c1833f35f8c244d983e51f0d5c49c` |
| `sim/aggregates/graph_detector_headline.json` | Committed run aggregate: graph-valued detectors on the full headline design. | 6153994 | `2c50a0bb7ea612e83db773fb920b250a509ce9b68c0fb0aa3bf85b52a338310b` |
| `sim/aggregates/graph_detector_low_coverage_stability.json` | Committed run aggregate: low-coverage observed-graph stability diagnostic. | 51183 | `51ec98d02be2c3cf201b3087e53c6ee055e5e53bd837c4586bf48a94fbefbffa` |
| `sim/aggregates/graph_detector_stratified_calibration_diagnosis.json` | Derived aggregate: diagnosis of graph-detector calibration transfer across strata. | 53863 | `412595fef26078ffe744f974b45b85c9356d333f3a5a147fe3471b7ad7df3c3c` |
| `sim/aggregates/graph_detector_stratified_eval_c025_partial.json` | Committed run aggregate: graph-detector stratified evaluation at coverage 0.25 (partial block). | 218118 | `701016549e8e2af0467995b56cf1b89e7376d8ee5707d8d05c8397d3f862a3d2` |
| `sim/aggregates/graph_detector_stratified_eval_c050_partial.json` | Committed run aggregate: graph-detector stratified evaluation at coverage 0.5 (partial block). | 217452 | `c4edf0cf1af0305f9dbfefc79f8c416074172f6458c96653c217f8dc9943f7ae` |
| `sim/aggregates/graph_detector_stratified_eval_c075_partial.json` | Committed run aggregate: graph-detector stratified evaluation at coverage 0.75 (partial block). | 218055 | `bde64c6dd8c69548ca3557518cb609b8162476bbbf442260d8cb7c949ab02689` |
| `sim/aggregates/graph_detector_stratified_eval_c100_partial.json` | Committed run aggregate: graph-detector stratified evaluation at coverage 1.0 (partial block). | 217641 | `e63669536a6e3471751b9ee8e257693d5240f8b79112da4ad01540a44ac98ecc` |
| `sim/aggregates/graph_detector_stratified_independent_null_cell00.json` | Committed run aggregate: graph-detector stratified independent null, grid cell 00. | 12094 | `6ee67556acc5acadee7de10c1898a1c1e1179a4233a7d2d4e782de7710f9ac5e` |
| `sim/aggregates/graph_detector_stratified_independent_null_cell01.json` | Committed run aggregate: graph-detector stratified independent null, grid cell 01. | 12094 | `b633235f3195d54d8415e46269e41a4c6adf82fc9f2e27f956d855902dac00c9` |
| `sim/aggregates/graph_detector_stratified_independent_null_cell02.json` | Committed run aggregate: graph-detector stratified independent null, grid cell 02. | 12120 | `275c5e06f2a088efe0f7f21b8e163dc86d70443fe59354b1fa12e51d99c92930` |
| `sim/aggregates/graph_detector_stratified_independent_null_cell03.json` | Committed run aggregate: graph-detector stratified independent null, grid cell 03. | 12031 | `279459e0170c43cbd4412b5e5718ccd70a6a7914bfc7e99a65758bfbbe43fb91` |
| `sim/aggregates/graph_detector_stratified_independent_null_cell04.json` | Committed run aggregate: graph-detector stratified independent null, grid cell 04. | 12027 | `ab8188bc173648b15147c12a370fc9427fe8aa43545bafc7b00e4bffe8381f77` |
| `sim/aggregates/graph_detector_stratified_independent_null_cell05.json` | Committed run aggregate: graph-detector stratified independent null, grid cell 05. | 12059 | `a03ce05baedeb91e750e39223155d911e2c68b54773bf13fdb4f30833163a104` |
| `sim/aggregates/graph_detector_stratified_independent_null_cell06.json` | Committed run aggregate: graph-detector stratified independent null, grid cell 06. | 12237 | `72231c35cce10915609a76f7b576f20daa7ce1c9e3854faca747d9e534af41f9` |
| `sim/aggregates/graph_detector_stratified_independent_null_cell07.json` | Committed run aggregate: graph-detector stratified independent null, grid cell 07. | 12228 | `9260c84ef5679193913f6b15c0ca238efeb1eee12dbb5d6b03d7e673be97ea29` |
| `sim/aggregates/graph_detector_stratified_independent_null_cell08.json` | Committed run aggregate: graph-detector stratified independent null, grid cell 08. | 12273 | `79a9001ddd45514f4964a0e81276a0677de23138604d9d18fe833846d4954037` |
| `sim/aggregates/graph_detector_stratified_independent_null_cell09.json` | Committed run aggregate: graph-detector stratified independent null, grid cell 09. | 12076 | `e45e4115c9c27c59b8f0800be271a408f851536e0738dee0ed0f5dd6be17bb62` |
| `sim/aggregates/graph_detector_stratified_independent_null_cell10.json` | Committed run aggregate: graph-detector stratified independent null, grid cell 10. | 12090 | `84fbce8e150dc3a20ba127b37e47a3e156899ffd43a005243e0ed20e4f705a7d` |
| `sim/aggregates/graph_detector_stratified_independent_null_cell11.json` | Committed run aggregate: graph-detector stratified independent null, grid cell 11. | 12113 | `ee6e3ede3cc2a085bbf92b5f42cd48037d15889a3d5621c0fb546889577d3a62` |
| `sim/aggregates/graph_detector_stratified_independent_null_cell12.json` | Committed run aggregate: graph-detector stratified independent null, grid cell 12. | 12027 | `0ab1930133e6ef9c93f218cc5b695ef55a2fb7a2b7a80141c9d1a523d0d44ccc` |
| `sim/aggregates/graph_detector_stratified_independent_null_cell13.json` | Committed run aggregate: graph-detector stratified independent null, grid cell 13. | 12029 | `e32b92d9545972fa72fd155cdd5f2f53213b19f2f9bf8c0e9877415d1eafb217` |
| `sim/aggregates/graph_detector_stratified_independent_null_cell14.json` | Committed run aggregate: graph-detector stratified independent null, grid cell 14. | 12058 | `37de7fa72e29acc5489fb37d77f2eca033dd428cc4d54ea6d71ee49a94b95552` |
| `sim/aggregates/graph_detector_stratified_independent_null_cell15.json` | Committed run aggregate: graph-detector stratified independent null, grid cell 15. | 12236 | `673b7b52fec9b5b08ea38df351be7058327981e5ff865b3ed71cfed810b31f82` |
| `sim/aggregates/graph_detector_stratified_independent_null_cell16.json` | Committed run aggregate: graph-detector stratified independent null, grid cell 16. | 12228 | `569e0f660c1fe7737fde46262b07fa306d4560dbd8914b320385147c12caa8dd` |
| `sim/aggregates/graph_detector_stratified_independent_null_cell17.json` | Committed run aggregate: graph-detector stratified independent null, grid cell 17. | 12250 | `07b49d85d81acb2385e7fc4d8b61b0dad9c755c14c77fad3afef84b3a6ecb514` |
| `sim/aggregates/graph_detector_stratified_independent_null_cell18.json` | Committed run aggregate: graph-detector stratified independent null, grid cell 18. | 12064 | `26c0d4bbf9725e34ef5c154d6c0c964a2a291d9ceb16de12d28bf1c9747d0a0e` |
| `sim/aggregates/graph_detector_stratified_independent_null_cell19.json` | Committed run aggregate: graph-detector stratified independent null, grid cell 19. | 12064 | `6ae0f6cd4f72e3a45f82760d46b8ca131e9afb5f690f59244c61d6b29567f41e` |
| `sim/aggregates/graph_detector_stratified_independent_null_cell20.json` | Committed run aggregate: graph-detector stratified independent null, grid cell 20. | 12092 | `b2b3049a0ae48f1e9f9d9813ad60ebcf264b01e2d9e93530e7ab2bec51f8909e` |
| `sim/aggregates/graph_detector_stratified_independent_null_cell21.json` | Committed run aggregate: graph-detector stratified independent null, grid cell 21. | 11999 | `1d40a1a0c131a11d9b83177f0ff469ef20854eea8548ec0c5e8d06044105dafe` |
| `sim/aggregates/graph_detector_stratified_independent_null_cell22.json` | Committed run aggregate: graph-detector stratified independent null, grid cell 22. | 12005 | `6719d1495bff749554d2cd8b33b8725231769d1623933af2572b35c0cef12223` |
| `sim/aggregates/graph_detector_stratified_independent_null_cell23.json` | Committed run aggregate: graph-detector stratified independent null, grid cell 23. | 12027 | `ffb93cf1589eac3d8c6af31dc9b71a7b874a10cb88d1ffc1e63ba8b021cd7cc7` |
| `sim/aggregates/graph_detector_stratified_independent_null_cell24.json` | Committed run aggregate: graph-detector stratified independent null, grid cell 24. | 12205 | `8cf70ef4218c904c618301e6e64806b0d944e2a00d7ab6365a5194627b5b74fc` |
| `sim/aggregates/graph_detector_stratified_independent_null_cell25.json` | Committed run aggregate: graph-detector stratified independent null, grid cell 25. | 12212 | `e51955bae9c5a217a90fdb7c1e282b8c9e025406081981ed3fb0d5bf996324a2` |
| `sim/aggregates/graph_detector_stratified_independent_null_cell26.json` | Committed run aggregate: graph-detector stratified independent null, grid cell 26. | 12242 | `f9724363258195e9190cc6d3160f902a29c4644f4d8f8b62592aa264bb05a4dd` |
| `sim/aggregates/graph_detector_stratified_independent_null_cell27.json` | Committed run aggregate: graph-detector stratified independent null, grid cell 27. | 12059 | `6f532e771ce0f186afe9ac2aef9ac9f1afc98a9c5bfcbef7b82ff87fe712296f` |
| `sim/aggregates/graph_detector_stratified_independent_null_cell28.json` | Committed run aggregate: graph-detector stratified independent null, grid cell 28. | 12047 | `2a48d27f26a84d27b4e673b8d28ba63922c8daa661ace21a13c07ffd7053e7bb` |
| `sim/aggregates/graph_detector_stratified_independent_null_cell29.json` | Committed run aggregate: graph-detector stratified independent null, grid cell 29. | 12091 | `6a4e4aa29c6a431d7b1314b30f3f50cef6e44e5d10a55d99135434fb71b4ff7a` |
| `sim/aggregates/graph_detector_stratified_independent_null_cell30.json` | Committed run aggregate: graph-detector stratified independent null, grid cell 30. | 11994 | `aaff8707bcf5c6beb39bd5913734c1868cc01d8274257fbbea6f8b847455af1e` |
| `sim/aggregates/graph_detector_stratified_independent_null_cell31.json` | Committed run aggregate: graph-detector stratified independent null, grid cell 31. | 12003 | `6baba210bdcf3c17d4f94797ca2e544fdfe7eeafeeebce4a0d671bf0c94b2b37` |
| `sim/aggregates/graph_detector_stratified_independent_null_cell32.json` | Committed run aggregate: graph-detector stratified independent null, grid cell 32. | 12029 | `f919b592c292626944fcb1f936b4281496aabdaafd9de9ae4a6b747bd0081275` |
| `sim/aggregates/graph_detector_stratified_independent_null_cell33.json` | Committed run aggregate: graph-detector stratified independent null, grid cell 33. | 12216 | `c0838c025a742f6e01e10cd2c95f3b716de5feaead7816e313feb1fbde1671f9` |
| `sim/aggregates/graph_detector_stratified_independent_null_cell34.json` | Committed run aggregate: graph-detector stratified independent null, grid cell 34. | 12212 | `15455b4a06c704d37b36fcb61d0225d842c50c9b7118c8447069dee9faa17d6d` |
| `sim/aggregates/graph_detector_stratified_independent_null_cell35.json` | Committed run aggregate: graph-detector stratified independent null, grid cell 35. | 12240 | `2e6de05d9abf3c40656bee206dc97e1a72793c6aac945157c18e75a76393ba97` |
| `sim/aggregates/graph_detector_stratified_independent_null_cell36.json` | Committed run aggregate: graph-detector stratified independent null, grid cell 36. | 12081 | `30c5dded4e285e31a4f984b40e04554a0db0066e49e5ed42553b1e4e7b15110a` |
| `sim/aggregates/graph_detector_stratified_independent_null_cell37.json` | Committed run aggregate: graph-detector stratified independent null, grid cell 37. | 12089 | `4d4fe6e38f65dfb1be27105fc10d1fcf94214e6774b49e3ae47d100c93487033` |
| `sim/aggregates/graph_detector_stratified_independent_null_cell38.json` | Committed run aggregate: graph-detector stratified independent null, grid cell 38. | 12121 | `767a4753e337795dc5dfa0d75084e10dd8b269d0ed095d5672905726cdf3eaf6` |
| `sim/aggregates/graph_detector_stratified_independent_null_cell39.json` | Committed run aggregate: graph-detector stratified independent null, grid cell 39. | 12034 | `2adc21abdc06a8d7c967ab06b81daf8032e7289689b9443f61c68abbe7506678` |
| `sim/aggregates/graph_detector_stratified_independent_null_cell40.json` | Committed run aggregate: graph-detector stratified independent null, grid cell 40. | 12032 | `58e93b9a096ae349aef342eae6c4f47775cc92caed8ee1638018f9ffbd75fb79` |
| `sim/aggregates/graph_detector_stratified_independent_null_cell41.json` | Committed run aggregate: graph-detector stratified independent null, grid cell 41. | 12057 | `a23af8185e3dd754f75911ea6d219b85805a16fad790315800ff5dbfdf7e9c57` |
| `sim/aggregates/graph_detector_stratified_independent_null_cell42.json` | Committed run aggregate: graph-detector stratified independent null, grid cell 42. | 12243 | `93b7b8ae080e46c8ceb61dbdf026e685a8ec2b7cde0a002c94360c63e6bb08a0` |
| `sim/aggregates/graph_detector_stratified_independent_null_cell43.json` | Committed run aggregate: graph-detector stratified independent null, grid cell 43. | 12241 | `3ea2cbe669ed2db011f2596e6a476aad237e729dad31dc53f4ad66fcce2c53c0` |
| `sim/aggregates/graph_detector_stratified_independent_null_cell44.json` | Committed run aggregate: graph-detector stratified independent null, grid cell 44. | 12265 | `cd3f841374e40890f278d1b84699403fb316d31d32d0f6648bb7d9b50f2ac857` |
| `sim/aggregates/graph_detector_stratified_independent_null_cell45.json` | Committed run aggregate: graph-detector stratified independent null, grid cell 45. | 12083 | `3e2cc77306ffd03d5bfbdee2514c891747c9eb3ba7af6e57719b459088b15711` |
| `sim/aggregates/graph_detector_stratified_independent_null_cell46.json` | Committed run aggregate: graph-detector stratified independent null, grid cell 46. | 12087 | `ea39cae00da2a2b1ded0f9f5a091e8b998918431aecebeaf837e43c744b7356c` |
| `sim/aggregates/graph_detector_stratified_independent_null_cell47.json` | Committed run aggregate: graph-detector stratified independent null, grid cell 47. | 12117 | `3a07706c39f0cedc34836589adb0213aa04a6c886abb3e470ca2906db81fe4d1` |
| `sim/aggregates/graph_detector_stratified_independent_null_cell48.json` | Committed run aggregate: graph-detector stratified independent null, grid cell 48. | 12023 | `dea9cc731d990d029ac7c756c4a292a090877a147f5b04a2feed57726757b3bd` |
| `sim/aggregates/graph_detector_stratified_independent_null_cell49.json` | Committed run aggregate: graph-detector stratified independent null, grid cell 49. | 12030 | `2a523bfdfe5b11342b6b227735379b3dfe6adc393a65bb7b302f4e33dd0c90e9` |
| `sim/aggregates/graph_detector_stratified_independent_null_cell50.json` | Committed run aggregate: graph-detector stratified independent null, grid cell 50. | 12051 | `65de11cc1d010d2c1c3b2dd7ef240f10d794e356cb135d4ca6609838bf53a62f` |
| `sim/aggregates/graph_detector_stratified_independent_null_cell51.json` | Committed run aggregate: graph-detector stratified independent null, grid cell 51. | 12237 | `da6001a85ea6cc15620c834aa684f30dff49505755b2c99f65337faddcc3d4ee` |
| `sim/aggregates/graph_detector_stratified_independent_null_cell52.json` | Committed run aggregate: graph-detector stratified independent null, grid cell 52. | 12237 | `48b441bd5663c855cd6ba266993e94ac67cfe10dd76a80159365d898c2d1a2a4` |
| `sim/aggregates/graph_detector_stratified_independent_null_cell53.json` | Committed run aggregate: graph-detector stratified independent null, grid cell 53. | 12277 | `7cbb9e8935b4c7061495c6958c2ebf81c49e30b6c164d38057919237cd886e41` |
| `sim/aggregates/graph_detector_stratified_independent_null_cell54.json` | Committed run aggregate: graph-detector stratified independent null, grid cell 54. | 12063 | `aa608799006e190e383266f53282a3c7b3b1e4ace020f3ab3636e510ab91beeb` |
| `sim/aggregates/graph_detector_stratified_independent_null_cell55.json` | Committed run aggregate: graph-detector stratified independent null, grid cell 55. | 12051 | `ddca4c47dfe86d4c2f4998786dfa0a231db681b856410f407722de8449707c4c` |
| `sim/aggregates/graph_detector_stratified_independent_null_cell56.json` | Committed run aggregate: graph-detector stratified independent null, grid cell 56. | 12086 | `8ed348675ff97d8afe6b231c7591963399d4225940ef85cde11bc3a1b6ab8606` |
| `sim/aggregates/graph_detector_stratified_independent_null_cell57.json` | Committed run aggregate: graph-detector stratified independent null, grid cell 57. | 11993 | `a5c8522a024adda410d39b33395fe3fd44bcdba2dead154cc7db847036ee6041` |
| `sim/aggregates/graph_detector_stratified_independent_null_cell58.json` | Committed run aggregate: graph-detector stratified independent null, grid cell 58. | 11991 | `40eef9e2fc150ab50fb3d041435a8968a9a85198330a33c4a45d76a3ae4c00de` |
| `sim/aggregates/graph_detector_stratified_independent_null_cell59.json` | Committed run aggregate: graph-detector stratified independent null, grid cell 59. | 12011 | `a76f5acaee061e6627dcff13465fd592952fe9670546c5e60b541ff774c7e8ff` |
| `sim/aggregates/graph_detector_stratified_independent_null_cell60.json` | Committed run aggregate: graph-detector stratified independent null, grid cell 60. | 12215 | `b3bbb999cba601bf8cbb469a6f8c7f25e1e6646eebec721bcecd52d9763cfa60` |
| `sim/aggregates/graph_detector_stratified_independent_null_cell61.json` | Committed run aggregate: graph-detector stratified independent null, grid cell 61. | 12202 | `cb3323cd3021e27c32a6d6ed317801eeee1304868255feb3b1166092d97cd511` |
| `sim/aggregates/graph_detector_stratified_independent_null_cell62.json` | Committed run aggregate: graph-detector stratified independent null, grid cell 62. | 12227 | `ee0b29360fde592daae4bf0e14948ed7ef53c8ea59812026080b060eab3dfeae` |
| `sim/aggregates/graph_detector_stratified_independent_null_cell63.json` | Committed run aggregate: graph-detector stratified independent null, grid cell 63. | 12053 | `89eb08f34ee269043f5cbb5c0219637ce7b2e2c329e45b26f22b2bebf93f0103` |
| `sim/aggregates/graph_detector_stratified_independent_null_cell64.json` | Committed run aggregate: graph-detector stratified independent null, grid cell 64. | 12058 | `6568d656f3c1cd71b0ea19e5460e82ae9738fe3269bcd8d9d3f5dcf04bbd2aa6` |
| `sim/aggregates/graph_detector_stratified_independent_null_cell65.json` | Committed run aggregate: graph-detector stratified independent null, grid cell 65. | 12082 | `ec6acd21a47816e3085bf12105d8c70f97e9af59a4f1f75f9b23247efd66623c` |
| `sim/aggregates/graph_detector_stratified_independent_null_cell66.json` | Committed run aggregate: graph-detector stratified independent null, grid cell 66. | 11993 | `8059abac7894d34afbc133ae930dbf0cc9f1b23dbbe72d0fc4aac40d3a3a5ad0` |
| `sim/aggregates/graph_detector_stratified_independent_null_cell67.json` | Committed run aggregate: graph-detector stratified independent null, grid cell 67. | 11999 | `8e38ff959350b7c5fb0fcb10b7c88af293615862f986d0080e4e0ba05293f176` |
| `sim/aggregates/graph_detector_stratified_independent_null_cell68.json` | Committed run aggregate: graph-detector stratified independent null, grid cell 68. | 12022 | `d906426ae4f1d7f859051858b11a5ff753e8a5cd8c75a373247b14cb9ca7cfd6` |
| `sim/aggregates/graph_detector_stratified_independent_null_cell69.json` | Committed run aggregate: graph-detector stratified independent null, grid cell 69. | 12201 | `04fda13745ea63c2b1a7d9fa2378ac4137af5e83fdfb9e7d61b17ea1ca1780dd` |
| `sim/aggregates/graph_detector_stratified_independent_null_cell70.json` | Committed run aggregate: graph-detector stratified independent null, grid cell 70. | 12198 | `29ce563ab18e3542b328c2c44bd0e05b35aa087c035308ce01586ce9f846e5d5` |
| `sim/aggregates/graph_detector_stratified_independent_null_cell71.json` | Committed run aggregate: graph-detector stratified independent null, grid cell 71. | 12240 | `b20a36a73751c2508853ff1c33c84c8d0b38b1a786680e59c8d18b59a4ae3e1e` |
| `sim/aggregates/graph_detector_stratified_null_c025_partial.json` | Committed run aggregate: graph-detector stratified null at coverage 0.25 (partial block). | 218262 | `152fd274b4d1d1855170c5bdd1916490eded946a21ac40ff3888f78adb0cf2ef` |
| `sim/aggregates/graph_detector_stratified_null_c050_partial.json` | Committed run aggregate: graph-detector stratified null at coverage 0.5 (partial block). | 217794 | `32189ad335af7e0de339da30bd0e6a0038e35d8ade25f587c317fb07ca8b58dc` |
| `sim/aggregates/graph_detector_stratified_null_c075_partial.json` | Committed run aggregate: graph-detector stratified null at coverage 0.75 (partial block). | 218100 | `0ac5a6d04716e4d7471260ddbeb2749c3f8d54f035eebf4c91f7a63fb6e345ff` |
| `sim/aggregates/graph_detector_stratified_null_c100_partial.json` | Committed run aggregate: graph-detector stratified null at coverage 1.0 (partial block). | 217605 | `e0e241b9c88ed35506dd37c738421cdff1c8dcc0e95e8c6cc6ae6cc5c04e78fb` |
| `sim/aggregates/missingness_arm_interval.json` | Derived aggregate: the MNAR-minus-MCAR achieved-edge difference and its 95% Wald interval at coverage 0.5, recomputed from ea_headline_recal per-run rows. | 2130 | `156914c51ac53fb067e2be85ef8a01c887b57f7c6406302f546c9e5d4af0c5fd` |
| `sim/aggregates/recal_analysis.json` | Derived aggregate: cross-campaign derivations from the recalibrated operating point. | 53202 | `e87a37f550223abd6eb79034c7c3a3fdafd2b6e66f084320a8b1838616c399cc` |
| `sim/aggregates/remediation_analysis.json` | Derived aggregate: level, shape and displacement summaries of the revised series. | 158203 | `435a86ed9391307b266d95d0a81d2102f5d182d965f4109ad714503b7e27258b` |
| `sim/aggregates/remediation_runs.json` | Committed run aggregate: raw remediation campaign runs. | 7878639 | `aaccdc54c2b2c82fa4c30b3c1b0e14eefa3970d17a693fed1cb7ef84bb1ad147` |
| `sim/aggregates/review_reanalysis.json` | Derived aggregate: THE REVISED INFERENTIAL RECORD produced by review_reanalysis.py — multiplicity family, Holm/BH survivors, recovery ratio and bias-reduction intervals. | 63308 | `00dee75b61acac630171fbbe3c9ef72e0d44e22f42b6c6113d3954d95b9d1c53` |
| `sim/aggregates/three_observable_block_a_calibration_20260903.json` | Committed run aggregate: block-A false-alarm matching for the three graph observables (2026-09-03). | 375304 | `40f14ee370c2567ec50c43893c941eb96c8f48bd25fda181984b43bbb6be728c` |
| `sim/aggregates/three_observable_block_b_results_20260903.json` | Committed run aggregate: block-B evaluation of the three graph observables (2026-09-03). | 3103004 | `0a45f6688eed4c55d42ff04b881893f34027c2376d8767315fb033a23eec3b29` |
| `sim/aggregates/tranche2_runs.json` | Committed run aggregate: follow-up tranche of knee-timing and reporting-noise variants. | 449263 | `30532ecd3a518fee806830af20cee305e5edc28e84ad2e0fbe1859d8beaea52c` |
| `sim/calibrate_class.py` | Both-class variant calibration. Writes configs/detection_rule_bothclass.json. | 8222 | `1bfbb490621a19429145fcef7b3877f401188f20a43cb323b0407623ef8f155c` |
| `sim/calibrate_recal.py` | Dedicated-null recalibration on 600 control runs (seed block 950000+). Writes configs/detection_rule_recal.json, the operating point the paper reports. | 13793 | `af61e22e2edb0eaae0cd6cede4192f05abc1b73f4460a281b1692cfd898d39df` |
| `sim/calibrate_rule.py` | Campaign-1 threshold calibration. Writes configs/detection_rule.json; never hand-edited. | 11113 | `9465448971bca5414fb04f6fa25d6f41230aaca8cb61ea1132415a6bad4c2a6f` |
| `sim/chi2.py` | Minimal chi-squared distribution support, implemented without scipy. | 3967 | `6f562f7470a70a5a1e2bfa55e8464eca8347d8a403bdb37645b158bcbb5b52d9` |
| `sim/claims.json` | The legacy claim binding: 128 historical ledger claims, each tied to the aggregate cell that produces it. | 27243 | `5d798bc455891c42347a9fa7ba1bcab1fda937c4108fd4432e1780a92567534d` |
| `sim/cluster_lead.py` | Cluster-level lead-time analysis over independent degradation groups. | 10801 | `a361f5493a3aac11b94b57761425691eee68788ca9086b6f80965f1de7220164` |
| `sim/cluster_matching.py` | Group construction and matching for the cluster analyses. | 6035 | `0355eea65241246d2f8f42f3f57b7350b8a4e5d7339db2ea703edf8bcaae1fe8` |
| `sim/cluster_recheck.py` | Independent recheck of the cluster construction. | 4818 | `aadae6f808a1d0e1fde14deb3351e56f276e3558a81f340592325902ed97cc8c` |
| `sim/configs/detection_rule.json` | Frozen campaign-1 detection rule, written by calibrate_rule.py. | 1354 | `3c8f61dd8e0adec755496ca0d538fcc74072e44888d1e22094901dc02795d2ba` |
| `sim/configs/detection_rule_bothclass.json` | Frozen both-class detection-rule variant (formation plus escalation latency). | 1874 | `0dce33516b19665d0bbbb2b8a5b27e0889b297c227192a2e662f050805820cc8` |
| `sim/configs/detection_rule_obs_giant_frac.json` | Frozen detection rule for the giant frac graph observable. | 1567 | `1d246da2d57c5121671142a7587e4190e59398677d1517be37bb02a479359e50` |
| `sim/configs/detection_rule_obs_giant_frac_split_20260903.json` | Frozen split-calibration detection rule for the giant frac graph observable (2026-09-03 block split). | 1228 | `a6b0e8790f91cba6991d6f3b3bb9939c0f7b739623e5c8baf4e60b362c35de89` |
| `sim/configs/detection_rule_obs_mean_degree.json` | Frozen detection rule for the mean degree graph observable. | 1572 | `ec4406c8db5d05f6a9f1bf90a79eb68c44dd0f2f6d9f0532788c2d4e0307e450` |
| `sim/configs/detection_rule_obs_mean_degree_split_20260903.json` | Frozen split-calibration detection rule for the mean degree graph observable (2026-09-03 block split). | 1228 | `cee73ec529dbfb9d527d68e6f99aee934f12040e5da1b012f8174bdcd0db079d` |
| `sim/configs/detection_rule_obs_susceptibility_split_20260903.json` | Frozen split-calibration detection rule for the susceptibility graph observable (2026-09-03 block split). | 1238 | `2f3211fc3c6659c41ccf669b2beeda902e1e88b2dba49f950b40d82f993f913f` |
| `sim/configs/detection_rule_recal.json` | Frozen recalibrated detection rule — the operating point the paper’s parameter table reports. | 2027 | `b40a6ed0b71dd137907269649fc64e254af7d92a215b62e9255c908c2dae7e8b` |
| `sim/emit_missingness_arm_interval.py` | Persists the MNAR-minus-MCAR achieved-edge interval the manuscript prints as an aggregate cell, recomputed from committed per-run rows. Changes no reported value. | 9181 | `87924eb85ca42632561a4fec5c14d465321b95e6488b27f64664827d191f2e67` |
| `sim/figures.py` | Figure generation for every published figure and its SVG source. | 28657 | `8dc28571be7c0a7896d864e16fa305718cfc0ffae4636df6ff3df3153001b1ae` |
| `sim/figures/fig1_ea_headline_mcar.png` | PNG render of the headline alarm-time gain under independent (MCAR) consent. | 126801 | `0760da3e04066e4bb7d749571268124849a7a4b95d258d5e98dfafecd707de9f` |
| `sim/figures/fig1_ea_headline_mcar.svg` | SVG source of the headline alarm-time gain under independent (MCAR) consent. | 67269 | `d17929c5fab2528a537442053fd5250b87f3a44bfa48af5b32130d564cbe399e` |
| `sim/figures/fig1_ea_headline_mcar_recal.png` | PNG render of the headline alarm-time gain under independent (MCAR) consent (recalibrated operating point). | 148053 | `a2eade4d7be4bd4144569146f69c23709e6cc038be059ecd1f45521528788f69` |
| `sim/figures/fig1_ea_headline_mcar_recal.svg` | SVG source of the headline alarm-time gain under independent (MCAR) consent (recalibrated operating point). | 79028 | `a8f8591f237d9846b736d8fb7375611b0d268d74b0714ee26ce21960bfcd3da2` |
| `sim/figures/fig1_ea_headline_mnar.png` | PNG render of the headline alarm-time gain under centrality-dependent (MNAR) consent. | 119655 | `ab8628aa01b9f28a7b05b5c2ba8b50fa5bddb28fcf770a4a859eb32970ebfe8a` |
| `sim/figures/fig1_ea_headline_mnar.svg` | SVG source of the headline alarm-time gain under centrality-dependent (MNAR) consent. | 66949 | `881ba183f1081eeeabab67b0a2413b87f0507479d7efee934498a82f756144b5` |
| `sim/figures/fig1_ea_headline_mnar_recal.png` | PNG render of the headline alarm-time gain under centrality-dependent (MNAR) consent (recalibrated operating point). | 149968 | `e603073f3d0e5722a105f3acf249fbf234f57809be71e8a1a06824e1026b900b` |
| `sim/figures/fig1_ea_headline_mnar_recal.svg` | SVG source of the headline alarm-time gain under centrality-dependent (MNAR) consent (recalibrated operating point). | 79993 | `49c1e54ad5896301b1c34627564ac79ce45af2dad78a904a4073fc5a582249b0` |
| `sim/figures/fig2_confusion_mnar_red.png` | PNG render of the alarm confusion series under centrality-dependent consent, correlated (red) error. | 125032 | `d1ca72b46fc6f5e2c6ad8405addf788edc55ee7c310a0a5bdf2536710b5cb9bc` |
| `sim/figures/fig2_confusion_mnar_red.svg` | SVG source of the alarm confusion series under centrality-dependent consent, correlated (red) error. | 70004 | `cb46cb18968b644d67e4e5cf04a73c2c87530a844bbab210355d391448655c49` |
| `sim/figures/fig2_confusion_mnar_red_recal.png` | PNG render of the alarm confusion series under centrality-dependent consent, correlated (red) error (recalibrated operating point). | 146491 | `c0d1b754698f216db572ba62b7790351c22d922488a55dbe66c332fb280e80d0` |
| `sim/figures/fig2_confusion_mnar_red_recal.svg` | SVG source of the alarm confusion series under centrality-dependent consent, correlated (red) error (recalibrated operating point). | 77319 | `152dd39f08aed70d85a1ecfc5c5c41708d9ac7795d0cc9f2cd7472e6318326b6` |
| `sim/figures/fig2_confusion_mnar_white.png` | PNG render of the alarm confusion series under centrality-dependent consent, white error. | 137592 | `b89d562af520afff059093005d44e0b2b4f16d30b4aae3a75ff4b7f514342d44` |
| `sim/figures/fig2_confusion_mnar_white.svg` | SVG source of the alarm confusion series under centrality-dependent consent, white error. | 72396 | `5a27e6aea35b9e02699e30bf5b136fbfc4c6ff6b9382cef6a96eafa1da39a2d9` |
| `sim/figures/fig2_confusion_mnar_white_recal.png` | PNG render of the alarm confusion series under centrality-dependent consent, white error (recalibrated operating point). | 161832 | `3ce8a924163db28ae180ea7414b716d9ee41dbc42ce9706815f2c0935b0b2f9a` |
| `sim/figures/fig2_confusion_mnar_white_recal.svg` | SVG source of the alarm confusion series under centrality-dependent consent, white error (recalibrated operating point). | 81571 | `3c711ea23eb84550d7aa024626f3da234c608e612e1c9a0b8bc9e58febee55b5` |
| `sim/figures/fig3_eb_joint.png` | PNG render of the joint coverage-by-coupling grid. | 142904 | `40eed04d6cf6af431359afa57278e9dc5ca5f907999d158ad50aeb027a1ec010` |
| `sim/figures/fig3_eb_joint.svg` | SVG source of the joint coverage-by-coupling grid. | 65179 | `f768fd5d6aa445442e993c97419b8c5fe97b9e943e96550827b065c63177749d` |
| `sim/figures/fig3_eb_joint_n500.png` | PNG render of the joint coverage-by-coupling grid at 500 matched pairs per cell. | 140140 | `a61c9ca7556db7506704a88423f01d90da5f897251fd27c17bf3fd781ff586e7` |
| `sim/figures/fig3_eb_joint_n500.svg` | SVG source of the joint coverage-by-coupling grid at 500 matched pairs per cell. | 64477 | `9b1138169ea00bea9d832cc02c1a58f7937f76313fd23af978b5260164e5e47f` |
| `sim/figures/fig4_isoline_recal.png` | PNG render of the retained coverage figure: alarm-time gain across coverage and noise severity, with the zero-gain iso-line (recalibrated operating point). | 129807 | `e31edb3d6017d938e4da6e773a9d53ac5bef7f2cfb736c91377521f78ed5947d` |
| `sim/figures/fig4_isoline_recal.svg` | SVG source of the retained coverage figure: alarm-time gain across coverage and noise severity, with the zero-gain iso-line (recalibrated operating point). | 116041 | `3067ae0833a45b5eddaae4bd53a34cb31d80db1532522db63f7c727042133c2a` |
| `sim/figures/fig5_calibration_null.png` | PNG render of the calibration-null false-alarm distribution. | 132535 | `4a5fdbb11d82c9c01b65870981262847726bc69d441eac966eb6a748ccd0265a` |
| `sim/figures/fig5_calibration_null.svg` | SVG source of the calibration-null false-alarm distribution. | 69306 | `6f7f20cd252a1fa8bccb8dd26866c321a85604ae4309288527cce7e47f33f571` |
| `sim/figures/fig6_eb_coupling_forest.png` | PNG render of the coupling-contrast forest plot. | 103630 | `410ede43ed42a23674e49c5ce2ee360e2c95894374d3a66a580d7d5b8ef3b088` |
| `sim/figures/fig6_eb_coupling_forest.svg` | SVG source of the coupling-contrast forest plot. | 68857 | `3a0d0a89274760e64793d3af21d5da7b71f3f69a78fb025a781b14f4f508229e` |
| `sim/figures/figA1_gate1.png` | PNG render of the bifurcation-validation appendix figure. | 151885 | `e7477cf0966bfaee7f84acadb2b5f3642ef0a987ba2f1c5aaf31f988d8059481` |
| `sim/figures/figA1_gate1.svg` | SVG source of the bifurcation-validation appendix figure. | 83944 | `42ad354123dcd81d3fc9da9294110dd5ff5c5b574814150fe883d212911a3765` |
| `sim/history/RESULTS-before-review-2026-09-11.md` | The pre-review results record, retained verbatim as provenance so the 2026-09-11 correction is inspectable. | 121711 | `ab28c722f7d2599f261a7b996dbb81e087c54f10b2d3bc87bf04bc9d6a09968e` |
| `sim/indicators.py` | The observation instrument: variance-family early-warning statistics, Kendall tau, AUC and the logistic consent-propensity fit, all without scipy. | 15008 | `bb723b00e0aa4ff634f63e18c9cf59581a9b83eb3e56b8e531d9400dd73cdddb` |
| `sim/model.py` | The agent model: team formation, repeat ties, attrition, the driver ramp, and the coverage-limited observation channel. | 46058 | `3e9e0aafa32e64e3855841b208e9d8b172007399c6cc9bd59ab63cf47e0ac5f1` |
| `sim/recal_analysis.py` | Cross-campaign derivations with no new runs: alarm-time decomposition, false-alarm intervals, seed-noise gauge, coupling interaction bound and anti-tracking. | 13948 | `c69dfc3c7bfc808547143f55b7b30f83efaa02e00627de8d002a1f942eb9afe0` |
| `sim/recompute_group_fpr.py` | Flat-control false-alarm rates by cluster bootstrap over 720 independent (coverage, arm, replicate) groups rather than over correlated runs. | 5384 | `c3ae96411a5ecc9c7c66cc004dc944836f79f75b414d89c692150cf0349e1da2` |
| `sim/remediation_analysis.py` | Analysis of the remediation runs: level, shape and displacement summaries of the revised series. | 40772 | `e7ef5795a3f352605ec8ac39b950ef3f619118035c8b5d6f200ec1570915d43a` |
| `sim/review_reanalysis.py` | THE REVISED INFERENTIAL RECORD. Enumerates the 64-comparison post-hoc family, applies Holm and Benjamini-Hochberg, and rebuilds the recovery-ratio and bias-reduction intervals on independent seed groups. Writes aggregates/review_reanalysis.json. | 11347 | `55819b786f7ee76974e562c27d5cbaf26b7b6248e8ba81b12c2eb31695560091` |
| `sim/run_aux_arms.py` | Auxiliary arms: revocation versus static-equivalent, and the propensity-correction arm, at fixed coverage and noise. | 22095 | `c0647a88fc1e4bbb71ef8e904f685142918fd97ce7343fcf67850a7d3ee95e7b` |
| `sim/run_coverage_finegrid.py` | Fine coverage grid (block B), the sixteen-cell grid used for coverage-heterogeneity tests. | 9989 | `ad710393ae8de1aaa9297a68665c69e9c8e20e72cfca5c79f9d9dce70dc5f52c` |
| `sim/run_ea_headline.py` | Headline coverage-and-noise campaign: coverage x noise colour x severity x consent arm, with paired controls. | 14079 | `192f1728f5e7fb26520d1e9b9c3414ce6b7e5d88a167e3e3548b20fa195a005b` |
| `sim/run_eb.py` | Coverage-and-coupling campaign: the joint coverage-by-coupling sweep with the coupling-severed control. | 16091 | `3694531a69aa0c3dadb10c31a1d93487198e5539c32fc2bbef5108aaf17301ae` |
| `sim/run_gate1_bifurcation.py` | Bifurcation validation: mean-field saddle-node, end-state bimodality, hysteresis and relaxation-time divergence. Must pass before any early-warning result is reported. | 20431 | `e2718e3b5058c95e1dc94f7e931e359133daadf67394b196856c7b891291f931` |
| `sim/run_graph_detector_campaign.py` | Graph-valued detector campaign: observed giant-component fraction, mean repeat-tie degree and percolation susceptibility, with calibration and stratified nulls. | 20821 | `3e3c2fec02052abfe99d32559628cd70851a3fd58cc9f24e87d6b4d20c4fa25d` |
| `sim/run_remediation.py` | Remediation campaign: the record-revision and bias-reduction arms under disclosed reduction and revocation. | 20332 | `521d177bf7271ac8e8088376194794907be27c2f66762196764c910feb6fa17e` |
| `sim/run_three_observable_campaign.py` | Three-observable campaign: block-A false-alarm matching and block-B evaluation for the three graph observables. | 13997 | `ecd3c1dceacd8f30605dc0f1465c076ccfb5605f98cc44e978ad18a3ffb82b4b` |
| `sim/run_tranche2.py` | Follow-up tranche: variant runs that report the headline result’s dependence on knee-timing and reporting-noise choices. | 32581 | `647d7594ccbbbf2c2d315a77036b813669f1cc874938ae5dedbf5ecc23216e84` |
| `sim/runner.py` | Execution harness. Runs a configuration and records the series a rule file is later scored against. | 3286 | `dfb4c267f37d86625e22294ea7bd4e9b40f0649f7f2f762d2a0262721f5991e7` |
| `sim/scoring.py` | The scoring rule: how a recorded run becomes an alarm time under a frozen detection rule. | 6180 | `a17397dfa515df07f5cb8c1c985d81bb8dbe6c9bd2dee94c6652aba35c7bec06` |
| `sim/svgwrite_min.py` | Dependency-free SVG writer, used as the figure backend when matplotlib is unavailable. | 7100 | `997fac035231bf8f004f85ec7d5f6fadf25e081e194920e9e047276ba7ac6ac0` |
| `sim/test_chi2.py` | Unit tests for chi2.py against known values. | 4430 | `2d16a4b89ffe94476f32b44ab25d572825501c9f32b6d169bd94afc694f1d7e7` |
| `sim/test_review_reanalysis.py` | Unit tests for the multiplicity correction and bootstrap constructions in review_reanalysis.py. | 1379 | `6b5ed9bc66ccb0649d8f2d49e2c66761e387f1479cc9003b0d35dd87f32aeef1` |
| `sim/verify_reported_numbers.py` | The legacy verifier: 3,873 checks re-deriving the earlier record from the committed aggregates. Exits non-zero on any mismatch. It binds the legacy record, not the revised manuscript. | 57480 | `fe0b515c6643dbd65912805f11662a5d7cfd6138e9dc69215f790bf091e6ac70` |
