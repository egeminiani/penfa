## Release summary

In this version (v0.1.1), I have fixed the vignette issue causing binary 
build failures after first CRAN release of penfa v0.1.0.


## Summary

penfa: a package for fitting single- and multiple-group penalized factor models 
via a trust-region algorithm with integrated automatic multiple tuning parameter 
selection.


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
checking CRAN incoming feasibility ... NOTE
Maintainer: 'Elena Geminiani <geminianielena@gmail.com>'

Days since last update: 3

```

* As requested by CRAN, I am submitting a new version fixing the binary build 
failures occurring after first package release to CRAN.
