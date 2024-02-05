if (!require("tidyverse", quietly = TRUE)) {
  install.packages("tidyverse")
}

if (!require("here", quietly = TRUE)) {
  install.packages("here")
}

render_reports <- function(meth_cut, analysis, tiling_window = "", tiling_step = ""){
  
  # Section 1: Create copies of main layout .qmd file
  # 1.a Create new names for .qmd copies by appending the species name to the
  # layout file name.
  file_in <- here::here(str_c("reports/report_lin_RRBS_F2_pancreas_meth_cut_", meth_cut, "_", analysis, tiling_window, "_", tiling_step, ".qmd"))
  
  # 1.b Create copies of the layout files using the modified file names
  file.copy(
    from = here::here("quarto/RRBS_lin_pancreas_xtro_male_F2.qmd"),
    to = file_in,
    overwrite = TRUE
  )
  
  # Section 2: Modify the copied .qmd file
  qmd_content <- readLines(file_in)
  for (i in seq_along(qmd_content)) {
    if (grepl("^  analysis:", qmd_content[i])) {
      qmd_content[i] <- paste0("  analysis: '", analysis, "'")
      break
    }
  }
  writeLines(qmd_content, file_in)
  
  # Section 3: Render reports using .qmd copies.
  quarto::quarto_render(
    input = file_in,
    execute_params = list(meth_cut = meth_cut,
                          analysis = analysis,
                          tiling_window = tiling_window,
                          tiling_step = tiling_step),
    execute_dir = here::here(),
    output_file = str_c("report_lin_RRBS_F2_pancreas_meth_cut_", meth_cut, "_", analysis, tiling_window, "_", tiling_step, ".html")
  )

}


# data frame with params
df <- data.frame(
  meth_cut = c(10, 15, 10, 20),
  analysis = c("CpGs", "CpGs", "tiles", "tiles"),
  tiling_window = c(1, 1, 100, 100),
  tiling_step = c(1, 1, 100, 100)
)

df <- data.frame(
  meth_cut = c(10, 15, 10, 20),
  analysis = c("CpGs", "CpGs", "tiles", "tiles"),
  tiling_window = c(1, 1, 100, 100),
  tiling_step = c(1, 1, 100, 100)
)

df <- data.frame(
  meth_cut = c(10, 15, 20),
  analysis = c("promoter", "promoter", "promoter"),
  tiling_window = c("", "", ""),
  tiling_step = c("", "", "")
)

# df <- data.frame(
#   meth_cut = c(10, 15),
#   analysis = c("tiles", "tiles"),
#   tiling_window = c(100, 100),
#   tiling_step = c(100, 100)
# )

# df <- data.frame(
#   meth_cut = c(10),
#   analysis = c("CpGs"),
#   tiling_window = c("", ""),
#   tiling_step = c("", "")
# )

for (i in 1:nrow(df)) {
  render_reports(meth_cut = df[i, "meth_cut"],
                 analysis = df[i, "analysis"],
                 tiling_window = df[i, "tiling_window"],
                 tiling_step = df[i, "tiling_step"])
  
  meth_cut <- df[i, "meth_cut"]
  analysis <- df[i, "analysis"]
  tiling_window <- df[i, "tiling_window"]
  tiling_step <- df[i, "tiling_step"]
  
  file.copy(from = here::here(str_c("report_lin_RRBS_F2_pancreas_meth_cut_", meth_cut, "_", analysis, tiling_window, "_", tiling_step, ".html")),
            to = here::here(str_c("reports/report_lin_RRBS_F2_pancreas_meth_cut_", meth_cut, "_", analysis, tiling_window, "_", tiling_step, ".html")))
  
  file.remove(from = here::here(str_c("report_lin_RRBS_F2_pancreas_meth_cut_", meth_cut, "_", analysis, tiling_window, "_", tiling_step, ".html")))
  }

