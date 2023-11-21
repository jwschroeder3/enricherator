# CHANGELOG

Significant updates the the Enricherator code base will be noted
in this file.

## [0.3.0] - dev

+ introduction of signal-to-noise allocation in a mixture model instead
of all reads having to come from the negative binomial process that fits
enrichments.

## [0.2.0] - 2023-09-18

### Fixed

+ v0.1.0 contained a bug that switched strand name
assignments when a minus strand file was listed first in the information
file. The bug did not affect our publications using Enricherator, because
we have always listed the plus stranded data first. The bug has been
fixed in v0.2.0.

## [0.1.0]

+ Release published with PMID: 37494439.

### Bug

+ See [version 0.2.0](#fixed) for note on bug involving strand assignments
in this version.
