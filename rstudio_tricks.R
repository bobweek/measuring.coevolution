# regexpr replace in rstudio


# to take an expression such as Power(AD,2) and replace with (AD^2) where A belongs to A-z and D belongs to 0-9
# use Power\(([A-z]\d),2\) as the search and ($1^2) as the replace.

# for Power(AD1,D2) where both D1 and D2 belong to 0-9, use Power\(([A-z]\d),(\d)\) and ($1^$2)