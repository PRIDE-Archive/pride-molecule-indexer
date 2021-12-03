spectra-indexer 
==============

This library index spectra lasrge scale in Elastic search for future visualization. The library contains a data model to store an spectra into an elastic document, also contains the information for retrival and query the elastic search.

## Data Model


The minimun infomarmation of an spectrum will be:

- usi: universal spectrum identifier : Please read the following [documentation](https://www.psidev.info/usi)
- masses: An array of peak masses
- intensities: An array of peak intensities
- charge: A charge state for the spectrum
- mz: mass to charge.


## Usage 

```maven
 <dependency>
     <groupId>io.github.bigbio.pgatk</groupId>
     <artifactId>spectra-indexer</artifactId>
     <version>${package.version}</version>
 </dependency>
 
```

 