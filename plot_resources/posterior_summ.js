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

CongraphPlotter.prototype.plot = function(cgraph, container, threshold_display) {
  var threshold = this._find_threshold(cgraph);
  d3.select(threshold_display).text(threshold.toFixed(2));
  var [nodes, edges] = this._make_graph(cgraph, threshold);
  this._draw(nodes, edges, container);
}

CongraphPlotter.prototype._find_threshold = function(cgraph) {
  var K = cgraph.length;
  var threshold = 1;

  for(var J = 1; J < K; J++) {
    var max_J = 0;
    for(var I = 0; I < K; I++) {
      if(cgraph[I][J] > max_J) {
        max_J = cgraph[I][J];
      }
    }
    if(threshold > max_J) {
      threshold = max_J;
    }
  }
  return threshold;
}

CongraphPlotter.prototype._make_graph = function(cgraph, threshold=0.2) {
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

CongraphPlotter.prototype._draw = function(nodes, edges, container) {
  var cy = cytoscape({
    container: document.querySelector(container),
    userZoomingEnabled: false,
    style: [{
      selector: 'node',
      style: {
        'background-color': '#428bca',
        'color': '#fff',
        'label': 'data(id)',
        'text-valign': 'center',
      },
    }, {
      selector: 'edge',
      style: {
        'width': 'data(thickness)',
        'opacity': 'data(weight)',
        'line-color': '#000',
        'target-arrow-color': '#000',
        'target-arrow-shape': 'triangle',
        // Necessary to allow arrowheads
        'curve-style': 'bezier',
        //'label': 'data(label)',
      },
    }],
  });

  nodes.forEach(node => {
    cy.add({
      group: 'nodes',
      data: { id: node.id },
    });
  });
  edges.forEach(edge => {
    const maxwt = 8;
    const minwt = 0.5;
    cy.add({
      group: 'edges',
      data: {
        id: edge.source + '_' + edge.target,
        source: edge.source,
        target: edge.target,
        weight: edge.weight,
        thickness: edge.weight*(maxwt - minwt) + minwt,
        label: edge.weight.toFixed(2),
      },
    });
  });

  cy.layout({name: 'fcose'}).run();
}
