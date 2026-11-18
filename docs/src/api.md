# API

```@meta
CurrentModule = Chopper
```

```@docs
Chopper
```

## Filtering files

```@docs
ChopperStats
chopper
```

## Filters and single reads

```@docs
ReadFilter
ReadStatus
read_segments
read_status
```

## Trimming strategies

```@docs
BestReadSegment
FixedCrop
NoTrim
SplitByLowQuality
TrimByQuality
TrimStrategy
trim_ranges
trim_ranges!
```

## Quality and GC content

```@docs
gc_content
mean_quality
phred_to_probability
```

## Contaminant screen

```@docs
ContaminantIndex
best_match_length
is_contaminant
```

## Modification tags

```@docs
ModTagError
rewrite_mod_tags
```

## Input

```@docs
open_fastq
```
