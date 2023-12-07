COMMENTS AND RESPONSES FROM Victoria Wimmer on DEC 6

If there are references describing the methods in your package, please
add these in the description field of your DESCRIPTION file in the form
authors (year) <doi:...>
authors (year) <arXiv:...>
authors (year, ISBN:...)
or if those are not available: <https:...>
with no space after 'doi:', 'arXiv:', 'https:' and angle brackets for
auto-linking.
(If you want to add a title as well please put it in quotes: "Title")

>>> Now adding a (In revisions) reference to the DESCRIPTION
>>> file in section Description

Please add \value to .Rd files regarding exported methods and explain
the functions results in the documentation. Please write about the
structure of the output (class) and also what the output means. (If a
function does not return a value, please document that too, e.g.
\value{No return value, called for side effects} or similar)
Missing Rd-tags:
      plot.dsem.Rd: \value
      predict.dsem.Rd: \value
      print.dsem.Rd: \value
      residuals.dsem.Rd: \value
      simulate.dsem.Rd: \value
      summary.dsem.Rd: \value
      vcov.dsem.Rd: \value

>>> We have added @return statements in the roxygen text that is
>>> compiled to build the .Rd files

Please always make sure to reset to user's options(), working directory
or par() after you changed it in examples and vignettes and demos. -->
nst/doc/vignette.R
e.g.:
oldpar <- par(mfrow = c(1,2))
...
par(oldpar)

>>> Thanks for pointing this out, it is now fixed in the one instances in the
>>> vignettes
