# Save GO Analysis Plots to Files

This function saves a list of ggplot2 objects to files in a chosen
format. The file names are constructed from the provided `list_name` and
the names of each plot in the list.

## Usage

``` r
save_plot(
  plot_list,
  list_name,
  path = "results/",
  device = "svg",
  height = 7,
  width = 15
)
```

## Arguments

- plot_list:

  A named list of ggplot2 plot objects.

- list_name:

  A character string to prefix the file names.

- path:

  A character string specifying the directory path where the plots will
  be saved. Default is `"results/"`.

- device:

  A character string specifying the output device/format. Default is
  `"svg"`. Other possible values include `"png"`, `"pdf"`, etc.

- height:

  A numeric value specifying the height of the saved plot (in inches).
  Default is 7.

- width:

  A numeric value specifying the width of the saved plot (in inches).
  Default is 15.

## Value

None. Files are written to disk.

## Details

The function iterates over each plot in `plot_list` and uses `ggsave` to
save each plot in the specified format.

## Examples

``` r
if (FALSE) { # \dontrun{
  # Assume up_plots is a list of ggplot objects from up_go()
  save_plot(up_plots, list_name = "up_GO", device = "png")
} # }
```
