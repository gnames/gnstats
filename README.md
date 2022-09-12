# GNstats

GNstats provides statistical taxonomic data for a collection of taxons.
It determines which Kingdom is the most prevalent for the collection,
provides Kingdom distribution, and finds the main taxon, which contains
at least 50% of all taxons.

## Usage

GNstats require objects that correspond to its `Hierarchy` interface. Such
interface has one method `Taxon()` which returns a slice of `Taxon` objects.
The `Taxon` object contains the following fields:

* ID string
* Name string
* RankStr string
* Rank Rank

The main constructor method calculates statistics on creation.

```go
hs := testData(t)
res := stats.New(hs, 0.7)
```
