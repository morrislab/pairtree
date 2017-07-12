function highlight_labels(cell, should_highlight) {
  cell = $(cell);
  var col_idx = cell.index();
  var row_idx = cell.closest('tr').index();
  // Don't highlight when hovering over first row or column, as these are labels.
  if(col_idx === 0 || row_idx == 0) { return; }

  var tbl = cell.closest('table');
  var row_label = tbl.find('tbody tr').eq(row_idx).find('td:first-child');
  var col_label = tbl.find('thead tr:first-child').children('th').eq(col_idx);

  var combined = row_label.add(col_label);
  combined.toggleClass('highlighted', should_highlight);
}

$(document).ready(function() {
  $('table.matrix td').hover(
    function() { highlight_labels(this, true);  },
    function() { highlight_labels(this, false); }
  );
});
