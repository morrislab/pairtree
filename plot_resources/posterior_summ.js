function PosteriorSumm() {
}

PosteriorSumm.prototype.plot = function(structs, container) {
  if(structs.length === 0) {
    return;
  }

  structs.forEach(function(struct, idx) {
    var row = d3.select(container).append('tr');
    var struct_id = 'structure_' + idx;

    row.append('td').text(idx);
    row.append('td').text((100 * struct.prob).toFixed(1) + '%');
    row.append('td').text(struct.nlglh.toFixed(3));
    row.append('td').text(struct.count);
    row.append('td').attr('id', struct_id);

    var root = 0;
    (new TreePlotter()).plot(root, struct.parents, struct.phi, struct.samples, '#' + struct_id);
  });
}

function CongraphPlotter() {
}

CongraphPlotter.prototype.plot = function(cgraph, container) {
  var [nodes, edges] = this._make_graph(cgraph);
  this._draw(nodes, edges, container);
}

CongraphPlotter.prototype._make_graph = function(cgraph, threshold=0.05) {
  var K = cgraph.length;
  var nodes = [];
  var edges = [];

  for(var I = 0; I < K; I++) {
    nodes.push({
      id: I,
      name: I,
      radius: 20,
    });
    for(var J = 0; J < K; J++) {
      if(cgraph[I][J] >= threshold) {
        edges.push({
          source: I,
          target: J,
          weight: cgraph[I][J],
        });
      }
    }
  }

  return [nodes, edges];
}

function _drag(simulation) {
  function dragstarted(d) {
    if (!d3.event.active) simulation.alphaTarget(0.3).restart();
    d.fx = d.x;
    d.fy = d.y;
  }

  function dragged(d) {
    d.fx = d3.event.x;
    d.fy = d3.event.y;
  }

  function dragended(d) {
    if (!d3.event.active) simulation.alphaTarget(0);
    d.fx = null;
    d.fy = null;
  }

  return d3.drag()
    .on('start', dragstarted)
    .on('drag', dragged)
    .on("end", dragended);
}

function _compute_coords(source, target) {
  const arrowhead_width = 14;
  const H = target.radius + arrowhead_width;
  const X = target.x - source.x;
  const Y = target.y - source.y;
  const theta = Math.atan(Math.abs(Y / X));

  var delta_x = H*Math.cos(theta);
  var delta_y = H*Math.sin(theta);
  if(X < 0) delta_x *= -1;
  if(Y < 0) delta_y *= -1;

  const coords = {
    x: target.x - delta_x,
    y: target.y - delta_y,
  };

  return coords;
}

CongraphPlotter.prototype._draw = function(nodes, edges, container) {
  var width = 800, height = 600;

  const simulation = d3.forceSimulation(nodes)
    .force('link', d3.forceLink(edges).id(d => d.id).distance(50))
    .force('charge', d3.forceManyBody().strength(-50))
    .force('center', d3.forceCenter(width / 2, height / 2));

  const svg = d3.select(container).append('svg')
    .attr('width', width)
    .attr('height', height);
  svg.append('svg:defs')
    .append('svg:marker')
    .attr('id', 'arrowhead')
    .attr('orient', 'auto')
    .attr('markerWidth', 10)
    .attr('markerHeight', 10)
    .attr('markerUnits', 'userSpaceOnUse')
    .attr('refX', 0)
    .attr('refY', 3)
    .append('svg:path')
    .attr('d', 'M0,0 L0,6 L9,3 z');

  const edge = svg.append('g')
    .attr('stroke', '#999')
    .selectAll('line')
    .data(edges)
    .join('line')
    .attr('stroke-opacity', d => d.weight)
    .attr('stroke-width', function(d) {
      var minwt = 0.5, maxwt = 4;
      return d.weight*(maxwt - minwt) + minwt;
    });

  edge.attr('marker-end', 'url(#arrowhead)');

  const node = svg.append('g')
    .attr('stroke', '#fff')
    .attr('stroke-width', 1.5)
    .selectAll('g')
    .data(nodes)
    .join('g');
  node.append('circle')
    .attr('r', d => d.radius)
    .attr('fill', '#ff0000')
    .call(_drag(simulation));
  node.append('text')
    .attr('font-size', '20')
    .attr('dominant-baseline', 'central')
    .attr('text-anchor', 'middle')
    .attr('fill', '#000000')
    .text(d => d.id);

  simulation.on('tick', () => {
    edge
      .attr('x1', d => d.source.x)
      .attr('y1', d => d.source.y)
      .attr('x2', d => _compute_coords(d.source, d.target).x)
      .attr('y2', d => _compute_coords(d.source, d.target).y);
    /*edge.attr('d', d3.linkHorizontal().source(function(d) {
      return [d.source.x + d.source.radius, d.source.y];
    }).target(function(d) {
      return [d.target.x - d.target.radius - arrowhead_width, d.target.y];
    }));*/

    node.attr('transform', d => 'translate(' + d.x + ',' + d.y + ')');
  });
}

/*CongraphPlotter.prototype._draw_old = function(nodes, edges, container) {
  var d3cola = cola.d3adaptor(d3).convergenceThreshold(0.1);
  var width = 960, height = 700;
  var outer = d3.select(container).append('svg')
    .attr('width',width)
    .attr('height',height)
    .attr('pointer-events','all');
  outer.append('rect')
    .attr('class','background')
    .attr('width','100%')
    .attr('height','100%')
    .call(d3.zoom().on('zoom', redraw));
  var vis = outer
    .append('g')
    .attr('transform', 'translate(250,250) scale(0.3)');

  function redraw() {
    vis.attr('transform', 'translate(' + d3.event.translate + ')' + ' scale(' + d3.event.scale + ')');
  }

  outer.append('svg:defs').append('svg:marker')
    .attr('id','end-arrow')
    .attr('viewBox','0 -5 10 10')
    .attr('refX',8)
    .attr('markerWidth',6)
    .attr('markerHeight',6)
    .attr('orient','auto')
    .append('svg:path')
    .attr('d','M0,-5L10,0L0,5L2,0')
    .attr('stroke-width','0px')
    .attr('fill','#000');

  d3cola.avoidOverlaps(true)
    .convergenceThreshold(1e-3)
    .flowLayout('x', 150)
    .size([width, height])
    .nodes(nodes)
    .links(edges)
    .jaccardLinkLengths(150);

  var link = vis.selectAll('.link')
    .data(edges)
    .enter().append('path')
    .attr('class', 'link');
  var margin = 10, pad = 12;
  var node = vis.selectAll('.node')
    .data(nodes)
    .enter().append('rect')
    .classed('node', true)
    .attr('rx',5)
    .attr('ry',5)
    .call(d3cola.drag);

  var label = vis.selectAll('.label')
    .data(nodes)
    .enter().append('text')
    .attr('class', 'label')
    .text(function (d) { return d.name; })
    .call(d3cola.drag)
    .each(function (d) {
      var b = this.getBBox();
      var extra = 2 * margin + 2 * pad;
      d.width = b.width + extra;
      d.height = b.height + extra;
    });

  var make_line = d3.line()
    .x(function (d) { return d.x; })
    .y(function (d) { return d.y; });

  var routeEdges = function () {
    d3cola.prepareEdgeRouting();
    link.attr('d', function (d) {
      return make_line(d3cola.routeEdge(d));
    });
  };

  d3cola.start(50, 100, 200).on('tick', function () {
    node.each(function (d) { d.innerBounds = d.bounds.inflate(-margin); })
    .attr('x', function (d) { return d.innerBounds.x; })
    .attr('y', function (d) { return d.innerBounds.y; })
    .attr('width', function (d) {
      return d.innerBounds.width();
    })
    .attr('height', function (d) { return d.innerBounds.height(); });

    link.attr('d', function (d) {
      var route = cola.makeEdgeBetween(d.source.innerBounds, d.target.innerBounds, 5);
      return make_line([route.sourceIntersection, route.arrowStart]);
    });

    label.attr('x', function (d) { return d.x })
      .attr('y', function (d) { return d.y + (margin + pad) / 2 });
  }).on('end', routeEdges);
}*/
