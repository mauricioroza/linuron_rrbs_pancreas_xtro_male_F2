if (!require("tidyverse", quietly = TRUE)) {
  install.packages("tidyverse")
}

if (!require("here", quietly = TRUE)) {
  install.packages("here")
}

render_reports <- function(meth_cut, group, analysis, tiling_window = "", tiling_step = ""){
  
  # Section 1: Create copies of main layout .qmd file
  # 1.a Create new names for .qmd copies by appending the species name to the
  # layout file name.
  file_in <- here::here(str_c("reports/report_IMZ_meth_cut_", meth_cut, "_", group, "_", analysis, tiling_window, "_", tiling_step, ".qmd"))
  
  # 1.b Create copies of the layout files using the modified file names
  file.copy(
    from = here::here("quarto/RRBS_IMZ_xtro_juv_liver_report.qmd"),
    to = file_in,
    overwrite = TRUE
  )
  
  # Section 2: Render reports using .qmd copies.
  quarto::quarto_render(
    input = file_in,
    execute_params = list(meth_cut = meth_cut,
                          group = group,
                          analysis = analysis,
                          tiling_window = tiling_window,
                          tiling_step = tiling_step),
    execute_dir = here::here(),
    output_file = str_c("report_IMZ_meth_cut_", meth_cut, "_", group, "_", analysis, tiling_window, "_", tiling_step, ".html")
  )

}


# data frame with params
df <- data.frame(
  meth_cut = c(10, 20, 10, 20),
  group = c("male", "male", "female", "female"),
  analysis = c("tiles", "tiles", "tiles", "tiles"),
  tiling_window = c(100, 100, 100, 100),
  tiling_step = c(100, 100, 100, 100)
)

df <- data.frame(
  meth_cut = c(10, 10),
  group = c("male", "female"),
  analysis = c("promoter", "promoter"),
  tiling_window = c("", ""),
  tiling_step = c("", "")
)

df <- data.frame(
  meth_cut = c(10),
  group = c("female"),
  analysis = c("CpGs"),
  tiling_window = c("", ""),
  tiling_step = c("", "")
)

for (i in 1:nrow(df)) {
  render_reports(meth_cut = df[i, "meth_cut"],
                 group = df[i, "group"],
                 analysis = df[i, "analysis"],
                 tiling_window = df[i, "tiling_window"],
                 tiling_step = df[i, "tiling_step"])
  
  meth_cut <- df[i, "meth_cut"]
  group <- df[i, "group"]
  analysis <- df[i, "analysis"]
  tiling_window <- df[i, "tiling_window"]
  tiling_step <- df[i, "tiling_step"]
  
  file.copy(from = here::here(str_c("report_IMZ_meth_cut_", meth_cut, "_", group, "_", analysis, tiling_window, "_", tiling_step, ".html")),
            to = here::here(str_c("reports/report_IMZ_meth_cut_", meth_cut, "_", group, "_", analysis, tiling_window, "_", tiling_step, ".html")))
  
  file.remove(from = here::here(str_c("report_IMZ_meth_cut_", meth_cut, "_", group, "_", analysis, tiling_window, "_", tiling_step, ".html")))
  }

