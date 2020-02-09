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
  const min_threshold = 0.01;
  const good_threshold = this._find_threshold(cgraph);

  var [nodes, edges] = this._make_graph(cgraph, min_threshold);
  this._draw(nodes, edges, container);
  this._config_threshold_chooser('#threshold_chooser', '#congraph_threshold', min_threshold, good_threshold);
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

CongraphPlotter.prototype._config_threshold_chooser = function(chooser, label_elem, min_thresh, good_thresh) {
  let self = this;
  const all_edges = this._cy.edges();

  let _update_thresh = function(new_thresh) {
    d3.select(label_elem).text(new_thresh);
    all_edges.restore();
    all_edges.filter('[weight <= ' + new_thresh + ']').remove();
    self._run_layout();
  };
  const rounded_thresh = Math.round(100*good_thresh)/100;
  _update_thresh(rounded_thresh);

  d3.select(chooser)
    .attr('min', min_thresh)
    .attr('max', 1.0)
    .attr('step', 0.01)
    .attr('value', rounded_thresh)
    .on('change', function(d) {
      const threshold = this.value;
      _update_thresh(threshold);
    });
}

CongraphPlotter.prototype._draw = function(nodes, edges, container) {
  let self = this;
  this._cy = cytoscape({
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
    self._cy.add({
      group: 'nodes',
      data: { id: node.id },
    });
  });
  edges.forEach(edge => {
    const maxwt = 8;
    const minwt = 0.5;
    self._cy.add({
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

  this._run_layout();
}

CongraphPlotter.prototype._run_layout = function() {
  this._cy.layout({name: 'fcose'}).run();
}
