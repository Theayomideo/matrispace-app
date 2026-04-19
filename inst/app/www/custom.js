/* MatriSpace - Custom JavaScript */

/* --- Initialization --- */

$(document).ready(function() {
  $(document).on('shiny:connected', function(event) {
    setTimeout(initializeCards, 100);
    initializeAllD3Visualizers();
  });

  if ($('#sel1').length > 0) { setTimeout(initializeCards, 100); }

  function initializeCards() {
    syncSel1Cards();
    syncSel2Cards();
  }

  $('#card_gene').on('click', function() { $('#sel1 input[value="matrisome gene"]').prop('checked', true).trigger('change'); });
  $('#card_signature').on('click', function() { $('#sel1 input[value="matrisome signature"]').prop('checked', true).trigger('change'); });
  $('#card_none').on('click', function() { $('#sel2 input[value="none"]').prop('checked', true).trigger('change'); });
  $('#card_gene2').on('click', function() { $('#sel2 input[value="matrisome gene"]').prop('checked', true).trigger('change'); });
  $('#card_any').on('click', function() { $('#sel2 input[value="any gene"]').prop('checked', true).trigger('change'); });
  $('#card_signature2').on('click', function() { $('#sel2 input[value="matrisome signature"]').prop('checked', true).trigger('change'); });

  function syncSel1Cards() {
    $('.feature-card:not(.secondary-card)').removeClass('selected');
    var sel1Value = $('#sel1 input:checked').val();
    if (sel1Value === 'matrisome gene') { $('#card_gene').addClass('selected'); }
    else if (sel1Value === 'matrisome signature') { $('#card_signature').addClass('selected'); }
  }
  function syncSel2Cards() {
    $('.secondary-card').removeClass('selected');
    var sel2Value = $('#sel2 input:checked').val();
    switch (sel2Value) {
      case 'none': $('#card_none').addClass('selected'); break;
      case 'matrisome gene': $('#card_gene2').addClass('selected'); break;
      case 'any gene': $('#card_any').addClass('selected'); break;
      case 'matrisome signature': $('#card_signature2').addClass('selected'); break;
    }
  }
  $('#sel1 input').change(syncSel1Cards);
  $('#sel2 input').change(syncSel2Cards);

  $(document).on('click', '.flip-card', function(event) { $(this).find('.flip-card-inner').toggleClass('is-flipped'); });

  $(document).on('click', '.accordion-button', function() {
    const fullText = $(this).text().trim();
    let panelName = 'unknown';
    if (fullText.includes('Load Data')) panelName = 'load_data';
    else if (fullText.includes('Profile Matrisome')) panelName = 'profile_matrisome';
    else if (fullText.includes('Select Feature')) panelName = 'feature_selection';
    else if (fullText.includes('Spatial Analysis')) panelName = 'spatial_analysis';
    else if (fullText.includes('Signature Analysis')) panelName = 'signature_analysis';
    Shiny.setInputValue('active_panel', panelName, { priority: "event" });
  });
});

/* --- D3 Visualization --- */

let vizInstances = {};

function getLayoutWidth(node) {
  if (!node) return 0;
  return node.clientWidth || node.offsetWidth || node.getBoundingClientRect().width || 0;
}

function getLayoutHeight(node) {
  if (!node) return 0;
  return node.clientHeight || node.offsetHeight || node.getBoundingClientRect().height || 0;
}

function initializeAllD3Visualizers() {
  $('.d3-viz-instance-container').each(function() {
    const prefix = $(this).data('prefix');
    if (!prefix) return;

    vizInstances[prefix] = {
      prefix: prefix,
      svg: d3.select(`#d3_svg_${prefix}`),
      tooltip: d3.select("body").append("div").attr("class", "tooltip").style("opacity", 0),
      cachedVizData: null,
      cachedImageDims: null,
      containerNode: this.querySelector(`#centerViz_${prefix}`),
      zoom: {
          panel: d3.select(`#zoom-panel_${prefix}`),
          container: d3.select(`#zoom-panel_${prefix} .zoom_container`),
          originalImg: d3.select(`#tissue_original_${prefix}`),
          vLine: d3.select(`#verticalLine_${prefix}`),
          hLine: d3.select(`#horizontalLine_${prefix}`),
      }
    };

    const observer = new ResizeObserver(entries => {
      if (entries[0].contentRect.width > 0) {
        drawVisualization(prefix);
      }
    });
    observer.observe(vizInstances[prefix].containerNode);
  });
}

function drawVisualization(prefix) {
  const instance = vizInstances[prefix];
  if (!instance || !instance.cachedVizData || !instance.svg) return;

  const { spots, imageUrl, clusterColors, imageDims, lowresSF, annotationName } = instance.cachedVizData;
  const svg = instance.svg;

  const containerWidth = getLayoutWidth(instance.containerNode);
  if (containerWidth === 0) return;

  const aspectRatio = imageDims.height / imageDims.width;
  const svgWidth = containerWidth;
  const svgHeight = containerWidth * aspectRatio;
  svg.attr("width", svgWidth).attr("height", svgHeight);

  svg.style("background", `url(${imageUrl}) no-repeat`).style("background-size", "100% 100%");
  svg.selectAll("circle").remove();

  const xScale = d3.scaleLinear().domain([0, imageDims.width]).range([0, svgWidth]);
  const yScale = d3.scaleLinear().domain([0, imageDims.height]).range([0, svgHeight]);

  const spotSize = Shiny.shinyapp.$inputValues[`spot_size_${prefix}`] || 3.5;
  const spotOpacity = Shiny.shinyapp.$inputValues[`spot_opacity_${prefix}`] || 0.8;
  const showBorders = Shiny.shinyapp.$inputValues[`spot_borders_${prefix}`] !== false;

  svg.selectAll("circle")
    .data(spots.filter(d => d.in_tissue == 1))
    .enter()
    .append("circle")
    .attr("class", "spot-circle")
    .attr("data-annotation", d => d[annotationName])
    .attr("cx", d => xScale(d.imagecol * lowresSF))
    .attr("cy", d => yScale(d.imagerow * lowresSF))
    .attr("r", spotSize)
    .style("fill", d => clusterColors[d[annotationName]] || "#CCCCCC")
    .style("opacity", spotOpacity)
    .style("stroke", showBorders ? "white" : "none")
    .style("stroke-width", showBorders ? 0.5 : 0)
    .on("mouseover", function(d) {
        instance.tooltip.style("opacity", 1);
        d3.select(this).style("stroke", "black").style("stroke-width", 2);
    })
    .on("mousemove", function(d) {
        instance.tooltip.html("<b>Annotation:</b> " + d[annotationName])
            .style("left", (d3.event.pageX + 15) + "px")
            .style("top", (d3.event.pageY - 30) + "px");
    })
    .on("mouseleave", function(d) {
        instance.tooltip.style("opacity", 0);
        d3.select(this).style("stroke", showBorders ? "white" : "none").style("stroke-width", showBorders ? 0.5 : 0);
    });

  if (prefix === "main") {
      const currentSelection = Shiny.shinyapp.$inputValues.selectedCellTypes_main;
      if (currentSelection) {
        updateSpotVisibility(prefix, currentSelection);
      }
  }

  setupZoom(prefix, imageUrl);
}

/* --- Zoom Panel --- */

function setupZoom(prefix, url) {
    const instance = vizInstances[prefix];
    const { svg, zoom } = instance;
    if (svg.empty() || zoom.panel.empty()) return;

    const svgNode = svg.node();
    const magnification = 2.0;

    const img = new Image();
    img.onload = () => { instance.cachedImageDims = { width: img.naturalWidth, height: img.naturalHeight }; };
    img.src = url;

    svg.on('mousemove.zoom', function() {
        if (!instance.cachedImageDims) return;

        zoom.originalImg.attr('src', url);
        zoom.container.classed('hidden', false);

        const svgWidth = parseFloat(svg.attr("width")) || getLayoutWidth(svgNode);
        const svgHeight = parseFloat(svg.attr("height")) || getLayoutHeight(svgNode);
        const zoomContainerSize = getLayoutWidth(zoom.panel.node());

        zoom.container.style("width", zoomContainerSize + "px").style("height", zoomContainerSize + "px");

        const originalWidth = instance.cachedImageDims.width;
        const originalHeight = instance.cachedImageDims.height;
        const zoomedImageWidth = originalWidth * magnification;
        const zoomedImageHeight = originalHeight * magnification;
        zoom.originalImg.style('width', zoomedImageWidth + 'px').style('height', zoomedImageHeight + 'px');

        const [mouseXInSvg, mouseYInSvg] = d3.mouse(this);
        const mouseXProportion = mouseXInSvg / svgWidth;
        const mouseYProportion = mouseYInSvg / svgHeight;
        const moX = -(mouseXProportion * zoomedImageWidth - zoomContainerSize / 2);
        const moY = -(mouseYProportion * zoomedImageHeight - zoomContainerSize / 2);
        zoom.originalImg.style('margin-left', moX + 'px').style('margin-top', moY + 'px');

        const centerX = zoomContainerSize / 2;
        const centerY = zoomContainerSize / 2;
        zoom.vLine.style('left', centerX + 'px');
        zoom.hLine.style('top', centerY + 'px');
    });

    svg.on('mouseout.zoom', function() {
        zoom.container.classed('hidden', true);
    });
}

/* --- Shiny Message Handlers --- */

Shiny.addCustomMessageHandler("clearViz", function(message) {
  Object.values(vizInstances).forEach(instance => {
      if (instance.svg) {
          instance.svg.selectAll("*").remove();
          instance.svg.style("background", "none");
      }
      instance.cachedVizData = null;
      instance.cachedImageDims = null;
  });
});

Shiny.addCustomMessageHandler("renderVizData", function(message) {
  const instance = vizInstances[message.target_id];
  if (!instance) return;
  instance.cachedVizData = message;
  drawVisualization(message.target_id);
});

function updateSpotVisibility(prefix, visibleAnnotations) {
    const instance = vizInstances[prefix];
    if (!instance || !instance.svg || !visibleAnnotations) return;
    instance.svg.selectAll(".spot-circle")
        .style("display", function() {
            const annotation = d3.select(this).attr("data-annotation");
            return visibleAnnotations.includes(annotation) ? "" : "none";
        });
}

Shiny.addCustomMessageHandler("updateVisibleAnnotations", function(message) {
    updateSpotVisibility(message.target_id, message.annotations);
});

Shiny.addCustomMessageHandler("updateSpotStyle", function(message) {
  const instance = vizInstances[message.target_id];
  if (!instance || !instance.svg) return;

  const { style, value } = message;
  const circles = instance.svg.selectAll(".spot-circle");

  switch (style) {
    case 'opacity':
      circles.style("opacity", value);
      break;
    case 'size':
      circles.attr("r", value);
      break;
    case 'borders':
      circles.style("stroke", value ? "white" : "none")
             .style("stroke-width", value ? 0.5 : 0);
      break;
  }
});
