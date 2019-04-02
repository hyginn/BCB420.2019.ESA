# `BCB420.2019.ESA`
# Notes on Pull Request Reviews

#### (**E**xploratory **S**ystems **A**nalysis **Tools for BCB420-2019**)
<!-- [![DOI](https://zenodo.org/badge/157482801.svg)](https://zenodo.org/badge/latestdoi/157482801) -->

&nbsp;


----

#### 1 General points

* Respect [**coding style**](http://steipe.biochemistry.utoronto.ca/abc/index.php/RPR-Coding_style)!
* Far too often I find pieces of code, or comments, or section names, or files that have nothing to do with the code, but were not removed from sample code that was copied and edited. **All of that needs to be fixed**. 
* Computations need to be **validated** for correctness at some point, in tests or in the vignette. Validation and testing are not the same. "Testing" checks whether you get the **expected behaviour** i.e. that there is not structural flaw in the algorithm's implementation. "Validation" confirms that the result is **correct** i.e. that there was no logical flaw in the algorithm's reasoning. 
* Headers need meaningful examples.
* Examples need comments on what they do.
* Examples with expensive computations or examples that need internet access should be placed in `\{dontrun ... }` blocks: e.g. 
```R
#' @examples
#' \dontrun{
#'myHGNC <- fetchData("HGNCreference")                   # fetch and assign data
#' }
```
* Mention the type and dimension of parameters and return values in the header. Type is one of `{numeric|integer|double|complex|logical|character|function}` ... or a mixed type such as `{data frame|list}`. The dimension is one of `{a|a vector|a matrix|a 3D matrix|...}`:
```R
#' @return (numeric)  A vector of the next six 6-of-49 lottery numbers or NULL.
```
* Be careful about return values: should they have the same length and order as the input? Eg. Should you return `NULL`, `character(0)`, `NA`, `NA_character`, `NaN`, or `""` for a no-match condition?
* Always `return()` explicitly, and always return something. For functions that are executed for their side-effects: `return(invisible(NULL))`.

&nbsp;

#### 2 Specific issues

**Issue**: _I am using this function or that data that was written by so-and-so_ ...

No need to mention a contributor's name separately since we are all collaborators in the same package. Authors can link to the specific function they are acknowledging instead. (e.g. `For details see \code{\link{<fetchData>}`). We can add some personalized attribution about ourselves into the Roxygen header of functions:
```R
#' @author <author name>, \email{<name>@@<host.domain>}
#' @seealso \code{\link{<name-of-related-function>}}
```
We can also include function authors' names in the LICENSE file.

----
**Issue**: Missing `@author` in the Roxygen header.

Use:
```R
#' @author <author name>, \email{<name>@@<host.domain>}
```


----



#### 3 References and Further Reading

(None)
&nbsp;

#### 4 Acknowledgements

(None)
&nbsp;

<!-- [END] -->
