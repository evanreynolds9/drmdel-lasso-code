# This script simply takes the simulation results and generates the text for the latex tables
source("R/init.R")

tol = 1e-8
basis_func = 12
latex_tables_file = file("latex_tables.txt", open = "w")
results_folder = "simulation_results"

header = c(
  "\\begin{table}[h]",
  "\\centering"
)

col_names = c(
  "\\begin{threeparttable}",
  "\\begin{tabular}{|l|r||r|r|r|r|r|}",
  "\\hline",
  "\\multirow{2}{*}{Estimator} & \\multirow{2}{*}{$n_k$} & \\multirow{2}{*}{Paths\\tnote{a}}",
  "& \\multicolumn{2}{c|}{$\\text{AIC}$} & \\multicolumn{2}{c|}{$\\text{BIC}$}  \\\\ \\cline{4-7}",
  "& & & Full\\tnote{a} & Sub\\tnote{a} & Full\\tnote{a} & Sub\\tnote{a} \\\\",
  "\\hline"
)

footer = c(
  "\\end{tabular}",
  "\\begin{tablenotes}",
  "\\small",
  "\\item[a] All figures in column are percentages.",
  "\\end{tablenotes}",
  "\\end{threeparttable}",
  "\\end{table}"
)

for(distribution in c("normal", "gamma", "lognormal")){
  if(distribution == "lognormal"){
    distribution_str = "log-normal"
  } else {
    distribution_str = distribution
  }
  
  caption = c(
    paste0("\\caption{Model selection results for ", distribution_str, " data}"),
    paste0("\\label{tab:", distribution, "_data}")
  )
  
  data_rows = c()
  
  for(est in c("reg", "adap", "aicbic")){
    
    if(est == "reg"){
      est_str = "\\multirow{5}{*}{$\\bm{\\hat{\\theta}}_{GL}$}"
    } else if(est == "adap"){
      est_str = "\\multirow{5}{*}{$\\bm{\\hat{\\theta}}_{AGL^*}$}"
    } else {
      est_str = "\\multirow{5}{*}{$\\bm{\\hat{\\theta}}_{MELE}$}"
    }
    
    data_rows = c(data_rows, est_str)
    
    for(n in c(250,500,1000,2500,5000)){
      sim_file = paste0(results_folder, "/", distribution, "_", n, "_", est, "_data.csv")
      if(est == "aicbic"){
        results = summariseAICBICSim(distribution = distribution, file_name = sim_file)
        run_entry = ""
      } else {
        results = summariseSim(distribution = distribution, file_name = sim_file, 
                               basis_func = basis_func, tol = tol)
        run_entry = format(round(getElement(results, "runs")*100, 1), nsmall = 1)
      }
      
      
      n_str = paste0(
        "& ", 
        paste(
          n, 
          run_entry,
          format(round(getElement(results, "AIC")*100, 1), nsmall = 1),
          format(round(getElement(results, "AIC_sub")*100, 1), nsmall = 1),
          format(round(getElement(results, "BIC")*100, 1), nsmall = 1),
          format(round(getElement(results, "BIC_sub")*100, 1), nsmall = 1),
          sep = " & "
        ),
        " \\\\"
      )
      
      data_rows = c(data_rows, n_str)
      
    }
    
    data_rows = c(data_rows, "\\hline")
    
  }
  
  # Write to the file
  writeLines(
    c(
      paste("Data for the", distribution, "distribution: "),
      "",
      header,
      caption,
      col_names,
      data_rows,
      footer,
      ""
    ), latex_tables_file
  )
  
}

close(latex_tables_file)