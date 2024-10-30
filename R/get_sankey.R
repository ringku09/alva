

enrich_sankey <- function(..., left_path = "KEGG", right_path = "RCTM",
                          font_size = 18, width = 1200, height=800,tag=NULL,
                          tag_size = 38, leg_size = 27,header_size = 40,
                          leg_text = "-log10(FDR)",
                          node_padding = 20, node_width = 30, header = NULL,
                          path_n = 10, organism = "rat", score_threshold = 200,
                          string_version = "11.5") {
  gene_set <- list2(...)
  if (length(gene_set) == 1 && is_bare_list(gene_set[[1]])) {
    gene_set <- gene_set[[1]]
  }

  if (is.null(names(gene_set))) {
    names(gene_set) <- paste("Set", LETTERS[1:length(gene_set)], sep = "_")
  }
  # gene_set <- lapply(probe_set, function(x) probes2genes(x))

  string_db <- setup_stringdb(organism = organism, score_threshold = score_threshold, version = string_version)

  gene_map <- gene_set # lapply(gene_set, function(x) string_db$map(x, "SYMBOL", removeUnmappedRows = TRUE))
  enrich_left <- lapply(gene_map, function(x) string_db$get_enrichment(x, category = left_path))
  enrich_right <- lapply(gene_map, function(x) string_db$get_enrichment(x, category = right_path))
  enrich_nleft <- lapply(enrich_left, function(x) if(nrow(x) > path_n) x[1:path_n,] else x)
  enrich_nright <- lapply(enrich_right, function(x) if(nrow(x) > path_n) x[1:path_n,] else x)
  left_all <- map2(lapply(enrich_nleft,"[", c("description","fdr")), names(enrich_nleft),
                    ~.x %>% mutate(geneset = .y, fdr = ceiling(-log10(fdr))))

  right_all <- map2(lapply(enrich_nright,"[", c("description","fdr")), names(enrich_nright),
                    ~.x %>% mutate(geneset = .y,fdr = ceiling(-log10(fdr))))
  left_df <- do.call(rbind, left_all)
  left_links <- left_df %>%
    relocate(source = description, target = geneset, value = fdr)
  right_df <- do.call(rbind, right_all)
  right_link <- right_df %>%
    relocate(source = geneset, target = description, value = fdr)
  links <- rbind(left_links, right_link) %>%
    mutate(group = as.factor(c(left_df$geneset, right_df$geneset)))
  n_set <- length(levels(links$group))
  group_let <- LETTERS[1:n_set]
  new_gr <- vector(length = length(links$group))
  for(i in 1: n_set) {
    new_gr[which(links$group == levels(links$group)[i])] = group_let[i]
  }
  nodes <- data.frame(
    name=c(as.character(links$source), as.character(links$target)) %>%
      unique()
  )
  links$group <- new_gr
  links$IDsource <- match(links$source, nodes$name)-1
  links$IDtarget <- match(links$target, nodes$name)-1 # "#FB9A99","#A6CEE3", "#B2DF8A", "#FDBF6F"
  nodes$group <- as.factor(c("my_nodes"))
  node_vec <- paste(group_let,collapse ="','")
  col_edg <- c("deepskyblue", "lightgreen", "#FB9A99","#CAB2D6", "#DDA0DD",
               "#CD853F", "#CAFF70", "#FFFF99", "#FDBF6F", "paleturquoise3")
  col_vec <- paste(c(col_edg[1:n_set], "gray"),collapse ="','")
  my_color <- htmlwidgets::JS(glue("d3.scaleOrdinal().domain(['{node_vec}']) .range(['{col_vec}'])"))
  leg_min <- width/2 - round(width*4.166/100)
  leg_max <- width/2 + round(width*4.166/100)
  leg_mid <- width/2
  x_min <- round(width*0.417/100)
  x_max <- round(width*0.833/100)
  p_min <- min(links$value)
  p_max <- max(links$value)
  if (is.null(header)) {
    header <- c(left_path, right_path)
  } else {
    header
  }
  sankey_lab <- list(header_size = paste0(header_size,"px"),
                     left_x = 0, right_x = width, left_header = header[1],
                     right_header = header[2])
  leg_data <- list(leg_min=leg_min,leg_max=leg_max, leg_mid = leg_mid, x_min=x_min,x_max=x_max,
                   p_min=p_min,p_max=p_max, leg_size = paste0(leg_size,"px"), leg_text = leg_text,
                   leg_wd = width*p_max/1000)
  midlab_data <- list(set_name = names(gene_set))

  p1 <- networkD3::sankeyNetwork(Links = links, Nodes = nodes, Source = "IDsource", Target = "IDtarget",
                     Value = "value", NodeID = "name",
                     colourScale = my_color, LinkGroup="group", NodeGroup="group",
                     fontSize = font_size ,width = width, height=height,
                     nodePadding = node_padding, nodeWidth = node_width, sinksRight = TRUE)
  # p2 <- htmlwidgets::onRender(p1, '
  # function(el,x, data) {
  #   //var cols_x = this.sankey.nodes().map(d => d.x).filter((v, i, a) => a.indexOf(v) === i).sort(function(a, b){return a - b});
  #   cols_x.forEach((d, i) => {
  #     d3.select(el).select("svg")
  #       .append("text")
  #       .attr("x", d+data.incm)
  #       .attr("y", 12)
  #       .attr("font-weight", "bold")
  #       .attr("text-anchor", [data.anchor[i]])
  #       .style("fill", "#696969")
  #       .style("font-size", data.header_size)
  #       .text([data.header[i]]);
  #   })
  # }', data = sankey_lab)
  p2 <- htmlwidgets::onRender(p1, '
  function(el,x, data) {
     var svg = d3.select("svg")
        svg.append("text")
        .attr("x", data.left_x)
        .attr("y", 9)
        .attr("font-weight", "bold")
        .attr("text-anchor", "start")
        .attr("alignment-baseline","central")
        .style("fill", "#696969")
        .style("font-size", data.header_size)
        .text(data.left_header);

       svg.append("text")
        .attr("x", data.right_x)
        .attr("y", 9)
        .attr("font-weight", "bold")
        .attr("text-anchor", "end")
        .attr("alignment-baseline","central")
        .style("fill", "#696969")
        .style("font-size", data.header_size)
        .text(data.right_header);
  }', data = sankey_lab)
  p3 <- htmlwidgets::onRender(p2,'
    function(el, x, data){
      var svg = d3.select("svg")
      // create data
      var dt = [{x: data.leg_min, y: data.p_min}, {x: data.leg_max, y: data.p_max}]

     // prepare a helper function
     var curveFunc = d3.area()
        .x(function(d) { return d.x })      // Position of both line breaks on the X axis
        .y1(function(d) { return d.y })     // Y position of top line breaks
        .y0(0)                            // Y position of bottom line breaks (200 = bottom of svg area)

     // Add the path using this helper function
     svg.append("path")
       .attr("d", curveFunc(dt))
       .attr("stroke", "#292929")
       .attr("fill", "#919191");
      svg.append("text").attr("x", data.leg_mid).attr("y", data.p_max).text(data.leg_text).style("font-size", data.leg_size).attr("alignment-baseline","hanging").attr("text-anchor", "middle")
      svg.append("text").attr("x", data.leg_min-data.x_max).attr("y", 0).text(data.p_min).style("font-size", data.leg_size).attr("alignment-baseline","hanging").attr("text-anchor", "end")
      svg.append("text").attr("x", data.leg_max+data.x_min).attr("y", 0).text(data.p_max).style("font-size", data.leg_size).attr("alignment-baseline","hanging")
   }',data=leg_data)

  p <- htmlwidgets::onRender(p3,'
       function(el, x, data){
          d3.select(el).selectAll(".node text").each((d, i, nodes) => {
            var currentNode = d3.select(nodes[i]);
            var nodeName = currentNode.text();
            if (data.set_name.includes(nodeName)) {
              nodes[i].setAttribute("x", x.options.nodeWidth/2);
              nodes[i].setAttribute("font-weight", "bold");
              nodes[i].setAttribute("text-anchor", "middle");}
     });
  }',data = midlab_data)

  # sn <- onRender(
  #   sn,
  #   '
  # function(el,x){
  # // select all our node text
  # d3.select(el)
  # .selectAll(".node text")
  # .filter(function(d) { return d.name.startsWith("Middle"); })
  # .attr("x", x.options.nodeWidth - 16)
  # .attr("text-anchor", "end");
  # }
  # '
  # )

  if(!is.null(tag)) {
    tag_level <- glue("font-family: Helvetica; font-size: {tag_size}px")
    p <- htmlwidgets::prependContent(p, htmltools::tags$h1(tag, style = tag_level))}
  return(p)
}




