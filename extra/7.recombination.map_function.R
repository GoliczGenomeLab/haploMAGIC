#!/usr/bin/env Rscript

#=================================================================================================
#title: 7.recombination.map_function.R
#description: plots recombination maps from the freq files containing the recombination frequencies calculated per SNP interval
#author: jmontero
#email: jose.a.montero-tena@ab.uni-giessen.de
#date: 2024-02-21
#version: 1.0.0
#notes: can be executed from 0.haplomagic.R
#       plots chromosome centromeric regions
#       plots SNP density as line marks
#=================================================================================================

PlotRecMap_Poisson = function (pop, chr, subset = 'RE') {
  # Read the centromere coordinates file
  centromereCoordinates = read.table('centromere_coordinates_rapeseed.txt', header = TRUE)
  
  # Extract Start(Mbp) and End(Mbp) values for the specified chromosome
  centromereStart = centromereCoordinates[as.integer(chr), 2]
  centromereEnd = centromereCoordinates[as.integer(chr), 3]
  
  # Read the frequency data
  outDF = read.table(paste0(pop, "_", chr, ".Poisson.freq"), colClasses = c(rep('character', 3), rep('numeric', 5)))
  
  # Find the top 3 maximum y-values and their corresponding x-values
  top_3_y_values <- head(sort(outDF[[ncol(outDF)]], decreasing = TRUE), 3)
  top_3_x_values <- outDF$pPos[order(outDF[[ncol(outDF)]], decreasing = TRUE)][1:3]/1e6
  
  # Creating a line plot with area under the curve filled in light blue
  plotRec = ggplot(data = outDF, aes(x = pPos/1000000, y = outDF[[ncol(outDF)]] )) +
    geom_area(fill = "#ADD8E6", alpha = 0.7) +  
    geom_line(color = "#00529B", linewidth = 0.5) +  
    geom_point(data = data.frame(x = outDF$pPos/1000000, y = min(outDF[[ncol(outDF)]])), aes(x, y), color = "black", size = 3, shape = 3) +  
    geom_point(data = data.frame(x = top_3_x_values, y = top_3_y_values), aes(x, y), color = "red", size = 3, shape = 19) +  
    geom_rect(aes(xmin = centromereStart, xmax = centromereEnd, ymin = -Inf, ymax = +Inf), fill = "grey", alpha = 0.01) + 
    geom_text(data = data.frame(x = top_3_x_values, y = top_3_y_values), aes(x, y, label = x), vjust = -0.5, color = "red", size = 3, angle = 90) +  
    labs(x = NULL, y = NULL,
         title = paste0("Recombination frequency (cM/Kbp) across chromosome ", chr, " (Mbp) in population ", pop, ".",
                        "\n", "Values transformed with Poisson correction.")) +
    scale_x_continuous(breaks = seq(0, round(max(outDF$pPos/1000000), 0), length.out = 11)) +  
    
    theme_minimal() +
    theme(
      panel.grid.major = element_line(color = "gray", linetype = "dotted"),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
      axis.line = element_line(color = "black", linewidth = 1),
      axis.text.x = element_text(angle = 0, hjust = 1),
      axis.title.y.right = element_blank(),
      axis.text.y.right = element_text(color = "black"),
      axis.ticks.y.right = element_line(color = "black"),
      axis.line.y.right = element_line(color = "black")
    ) +
    scale_y_continuous(position = "right")  
  
  # Save the plot as a PDF file
  print(plotRec)
  ggsave(paste0(pop, "_", chr, "_", subset, "_Poisson_RF.pdf"), plotRec, width = 8, height = 6)
}
PlotRecMap = function (pop, chr, subset, method) {
  # Read the centromere coordinates file
  centromereCoordinates = read.table('centromere_coordinates_rapeseed.txt', header = TRUE)

  # Extract Start(Mbp) and End(Mbp) values for the specified chromosome
  centromereStart = centromereCoordinates[chr, 2]
  centromereEnd = centromereCoordinates[chr, 3]

  #  outDF = read.table(paste0(pop, "_", chr, ".", method, ".freq"), colClasses = c(rep("character", 3), rep("numeric", 5)))
  outDF = fread(paste0(pop, "_", chr, ".", method, ".freq"), colClasses = c(rep('character', 3), rep('numeric', 5)))
  # Find the top 3 maximum y-values and their corresponding x-values
  top_3_y_values <- head(sort(outDF[[paste0("RE.RecFreq.cM.Kbp.", method)]], decreasing = TRUE), 3)
  top_3_x_values <- outDF$pPos[order(outDF[[paste0("RE.RecFreq.cM.Kbp.", method)]], decreasing = TRUE)][1:3]/1e6
  # Creating a line plot with area under the curve filled in light blue
  plotRec = 
    ggplot(data = outDF, aes(x = pPos/1000000, y = outDF[[paste0("RE.RecFreq.cM.Kbp.", method)]] )) +
    geom_area(fill = "#ADD8E6", alpha = 0.7) +  # Filling area under the curve in light blue with transparency
    geom_line(color = "#00529B", linewidth = 0.5) +  # Halved line width with a darker blue color
    geom_point(data = data.frame(x = outDF$pPos/1000000, y = min(outDF[[ncol(outDF)]])), aes(x, y), color = "black", size = 3, shape = 3) +  
    geom_point(data = data.frame(x = top_3_x_values, y = top_3_y_values), aes(x, y), color = "red", size = 3, shape = 19) +  # Mark top 3 y-values with red dots
    geom_rect(aes(xmin = centromereStart, xmax = centromereEnd, ymin = -Inf, ymax = +Inf), fill = "grey", alpha = 0.01) + 
    geom_text(data = data.frame(x = top_3_x_values, y = top_3_y_values), aes(x, y, label = x), vjust = -0.5, color = "red", size = 3, angle = 90) +  # Display y-values rotated 90 degrees
    labs(x = NULL, y = NULL,  # Blank axis titles
         title = paste0("Recombination frequency (cM/Kbp) across chromosome ", chr, " (Mbp) in population ", pop,
                        "\n", "METHOD: ", method, ". ")) +
    scale_x_continuous(breaks = seq(0, round(max(outDF$pPos/1000000), 0), length.out = 11)) +  # 10 ticks from 0 to max_x
    
    # Customize axis lines, background grids, and y-axis labels rotation
    theme_minimal() +
    theme(
      panel.grid.major = element_line(color = "gray", linetype = "dotted"),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
      axis.line = element_line(color = "black", linewidth = 1),
      axis.text.x = element_text(angle = 0, hjust = 1),  # Rotate x-axis labels
      axis.title.y.right = element_blank(),  # Blank y-axis label on the right side
      axis.text.y.right = element_text(color = "black"),  # Color y-axis text on the right side
      axis.ticks.y.right = element_line(color = "black"),  # Color y-axis ticks on the right side
      axis.line.y.right = element_line(color = "black")  # Color y-axis line on the right side
    ) +
    
    scale_y_continuous(position = "right")  # Move y-axis to the right side
  # Save the plot as a PDF file
  suppressWarnings(print(plotRec))
  ggsave(paste0(pop, "_", chr, "_", subset, "_", method, "_RF.pdf"), plotRec, width = 8, height = 6)  # You can adjust width and height as needed
}
