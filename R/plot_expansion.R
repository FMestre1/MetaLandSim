plot_expansion <-
function(expansion_object)
  {

  p1 <- gvisLineChart(expansion_object, 
                         xvar = "DISTANCE", 
                         yvar = "PROPORTION", 
                         options = list(title = "Dispersal Simulation", 
                                        width = 500, 
                                        height = 300, 
                                        curveType = "function", 
                                        legend = "none", 
                                        titleTextStyle = "{colour:'black', fontName:'Courier', fontSize:16}", 
                                        vAxis = "{title: 'proportion'}", 
                                        hAxis = "{title: 'distance (meters)'}", 
                                        series = "[{color: '#006400'}]", 
                                        backgroundColor = "#D1EEEE")
  )
  
  p1$html$caption <- paste("<div><span>Range expansion.</span><br />", 
                              sep = "")
  
  p1$html$footer <- paste("\n<!-- htmlFooter -->\n<span> \n  ",
                             R.Version()$version.string,
                             " &#8226; <a href=\"http://code.google.com/p/google-motion-charts-with-r/\">googleVis-", 
                             packageVersion("googleVis"),
                             "</a>\n  &#8226; MetaLandSim-",
                             packageVersion("MetaLandSim"),
                             "\n  &#8226; <a href=\"https://developers.google.com/terms/\">Google Terms of Use</a> &#8226; <a href=\"https://google-developers.appspot.com/chart/interactive/docs/gallery/linechart.html#Data_Policy\">Data Policy</a>\n</span></div>\n</body>\n</html>\n", 
                             sep="")

    
  plot(p1)

  }
