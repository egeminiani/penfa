## Resubmission

This is a resubmission. I have:

* Expanded the documentation to include all missing Rd-tags, explained 
  the function results and the structure of the output. 

* Removed the global `options(width=100)` in vignettes and only added the
  option temporarily to specific R chunks via `R.options = list(width = 100)`.


## Summary

First release of penfa, a package for fitting single- and multiple-group 
penalized factor models via a trust-region algorithm with integrated automatic 
multiple tuning parameter selection.


## Test environments

* local OS Windows 10 install, R 4.1.0
* Github Actions:
  - windows-latest (release)
  - macOS-latest (release)
  - ubuntu-20.04 (release)
  - ubuntu-20.04 (devel)
* win-builder (release, devel, oldrelease)


## R CMD check results

There were no ERRORs or WARNINGs. 

There was 1 NOTE on win-builder:

```
* checking CRAN incoming feasibility ... NOTE
Maintainer: 'Elena Geminiani <geminianielena@gmail.com>'

New submission

Possibly mis-spelled words in DESCRIPTION:
  Geminiani (23:26)
  al (23:39)
  et (23:36)
  mcp (24:62)
```

The above words in the DESCRIPTION are false positives, they are spelled 
correctly.

* This is a new release.
