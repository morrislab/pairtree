function PosteriorSumm() {
}

PosteriorSumm.prototype.plot = function(results, container) {
  if(results.structs.length === 0) {
    return;
  }

  results.structs.forEach(function(struct, idx) {
    var row = d3.select(container).append('tr');
    var struct_id = 'structure_' + idx;
    var tree_container = '#' + struct_id;

    row.append('td').text(idx);
    row.append('td').text((100 * struct.prob).toFixed(1) + '%');
    row.append('td').text(struct.nlglh.toFixed(3));
    row.append('td').text(struct.count);
    row.append('td').attr('id', struct_id);

    var root = 0;
    (new TreePlotter()).plot(
      root,
      struct.parents,
      struct.phi,
      struct.samples,
      results.samp_colours,
      null,
      false,
      tree_container,
    );
    // Resize tree to fit in table.
    d3.select(tree_container)
      .select('svg')
      .attr('width', '100%');
  });
}

function CongraphPlotter() {
  this._layout_options = {
    klay: {
      name: 'klay',
      animate: true,
      klay: {
        direction: 'RIGHT',
      },
    },

    fcose: {
      name: 'fcose',
    },
  };

  this._default_layout = this._layout = 'klay';
}

CongraphPlotter.prototype.plot = function(cgraph, container, threshold_display) {
  const self = this;
  const min_threshold = 0.01;
  const default_threshold = 0.05;

  const spanning_threshold = this._find_min_spanning_threshold(cgraph);
  if(spanning_threshold < min_threshold) {
    throw "spanning_threshold < min_threshold, implying some parts of the graph will always be disconnected";
  }
  d3.select('#spanning_threshold').text(Math.round(100*spanning_threshold) + '%');

  var [nodes, edges] = this._make_graph(cgraph, min_threshold);
  this._draw(nodes, edges, container);

  this._cy.ready(() => {
    self._config_layout_chooser('#layout_chooser');
    self._config_edge_weight_display();
    self._config_threshold_chooser('#threshold_chooser', '#congraph_threshold', min_threshold, default_threshold);
    self._config_exporters('#export_svg', '#export_png');
  });
}

CongraphPlotter.prototype._find_min_spanning_threshold = function(cgraph) {
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

CongraphPlotter.prototype._config_threshold_chooser = function(chooser, label_elem, min_thresh, default_thresh) {
  let self = this;
  const all_edges = this._cy.edges();

  let _update_thresh = function(new_thresh) {
    d3.select(label_elem).text(Math.round(100*new_thresh) + '%');
    all_edges.restore();
    const bad = all_edges.filter(E => E.data().weight < new_thresh);
    bad.remove();
    self._run_layout();
  };
  _update_thresh(default_thresh);

  d3.select(chooser)
    .attr('min', min_thresh)
    .attr('max', 1.0)
    .attr('step', 0.01)
    .attr('value', default_thresh)
    .on('change', function(d) {
      const threshold = this.value;
      _update_thresh(threshold);
    });
}

CongraphPlotter.prototype._config_layout_chooser = function(chooser) {
  var self = this;

  chooser = d3.select(chooser);
  for(const layout of Object.keys(this._layout_options).sort()) {
    const option = chooser.insert('option').text(layout);
    if(layout == this._default_layout) {
      option.attr('selected', 'selected');
    }
  }

  chooser.on('change', function(d) {
    const new_layout = this.value;
    self._layout = new_layout;
    self._run_layout();
  });
}

CongraphPlotter.prototype._config_exporters = function(svg_sel, png_sel) {
  let self = this;
  d3.select(svg_sel).on('click', function(d) {
    let svg = self._cy.svg();
    let blob = new Blob([svg], {type:"image/svg+xml;charset=utf-8"})
    saveAs(blob, 'congraph.svg');
  });
  d3.select(png_sel).on('click', function(d) {
    let blob = self._cy.png({
      output: 'blob',
      scale: 5,
    });
    saveAs(blob, 'congraph.png');
  });
}

CongraphPlotter.prototype._config_edge_weight_display = function() {
  // See https://stackoverflow.com/a/54556015
  // and https://cytoscape.org/cytoscape.js-popper/demo-tippy.html
  this._cy.edges().forEach(function(edge) {
    var ref = edge.popperRef();
    const dummy = document.createElement('div');

    edge.tip = tippy(dummy, {
      onCreate: function(inst) {
        inst.popperInstance.reference = ref;
      },
      lazy: false,
      trigger: 'manual',
      content: () => {
        const content = document.createElement('div');
        content.innerHTML = Math.round(100*edge.data().weight) + '%';
        return content;
      },
    });
  });

  this._cy.edges().bind('mouseover', (evt) => evt.target.tip.show());
  this._cy.edges().bind('mouseout',  (evt) => evt.target.tip.hide());
}

CongraphPlotter.prototype._draw = function(nodes, edges, container) {
  let self = this;
  this._cy = cytoscape({
    container: document.querySelector(container),
    wheelSensitivity: 0.1,
    userZoomingEnabled: true,
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
        'opacity': 'data(opacity)',
        'line-color': '#000',
        'target-arrow-color': '#000',
        'target-arrow-shape': 'triangle',
        // Necessary to allow arrowheads
        'curve-style': 'bezier',
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
    const max_opacity = 1.0;
    const min_opacity = 0.1;

    self._cy.add({
      group: 'edges',
      data: {
        id: edge.source + '_' + edge.target,
        source: edge.source,
        target: edge.target,
        weight: edge.weight,
        thickness: edge.weight*(maxwt - minwt) + minwt,
        opacity: edge.weight*(max_opacity - min_opacity) + min_opacity,
      },
    });
  });

  this._run_layout();
}

CongraphPlotter.prototype._run_layout = function() {
  this._cy.layout(this._layout_options[this._layout]).run();
}
