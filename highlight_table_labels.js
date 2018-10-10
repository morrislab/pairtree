function highlight_labels(cell, should_highlight) {
  cell = $(cell);
  var col_idx = cell.index();
  var row_idx = cell.closest('tr').index();

  var tbl = cell.closest('table');
  var row_label = tbl.find('tbody tr').eq(row_idx).find(':first-child');
  var col_label = tbl.find('thead tr:first-child').children().eq(col_idx);

  var combined = row_label.add(col_label);
  combined.toggleClass('highlighted', should_highlight);
}

$(document).ready(function() {
  $('table.matrix td').hover(
    function() { highlight_labels(this, true);  },
    function() { highlight_labels(this, false); }
  );
});
