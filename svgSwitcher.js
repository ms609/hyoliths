// Enable consensus tree interactions:
function showNodeNotes (e) {
  var target = e.target; // Which node triggered the click
  if (target.nodeName === 'tspan')
    $target = $(target);
    gp = $target.parent().parent();
    if (gp[0].nodeName === 'g' && gp.hasClass('node-g')) {
      nodeNumber = gp.data('node');
      $("#conTree g").removeClass("bold");
      gp.addClass("bold");
      document.getElementById('nodeDetails').style='display: none';
      $('div .nddtl').css('display', 'none');
      document.getElementById('nodeDetail' + nodeNumber).style='display: block';
  }
  return false;
}
$("#conTree").click(showNodeNotes)


// Hide details if not desired
function toggleSVG() {
  $('svg.tree').toggle();
  $('div.switcher').toggle();
  if (document.getElementById('showOrHide').innerHTML == 'hide') {
    document.getElementById('showOrHide').innerHTML = 'show'
  } else {
    document.getElementById('showOrHide').innerHTML = 'hide'
  }
}

function toggleDetails() {
  $this = $(this);
  index = $this.data('index');
  $('#stateNotes' + index).slideToggle();
  if ($this.html() == "[Hide details]") {
    // Leave tree, or it's too bouncy.
    $this.html("[Show details]");
  } else {
    $('svg.tree').eq(index).slideDown();
    $('div.switcher').eq(index).slideDown();
    $this.html("[Hide details]");
  }
}

// Prepare for character-by-character interactions
var latestKnownScrollY = 0,
  latestKnownY = 0,
	ticking = false;

function onScreen(e, viewTop, viewHeight) {
  var viewBottom = viewTop + viewHeight;
  var elemTop = $(e).offset().top;
  return ((elemTop <= viewBottom + 1000) && (elemTop >= viewTop - 1000));
}

function highlightTaxon (e) {
  var target = e.target; // Which node triggered the click
  $target = $(target);
  if (target.nodeName === 'text' && $target.hasClass('taxonLabel')) {
    taxonText = $(target).text();
    $withSameText = $('#reconstructions svg text:contains("' + taxonText + '")')
    if (localStorage.getItem(taxonText)) {
      localStorage.removeItem(taxonText)
      $withSameText.removeClass("highlightTaxon");
    } else {
      localStorage.setItem(taxonText, 1)
      $withSameText.addClass("highlightTaxon");
    }
  }
  return false;
}

function applyHighlighting () {
  if(localStorage.getItem($(this).text())) {
    $(this).addClass("highlightTaxon");
  }
}

function onScroll() {
  //Optimised via https://www.html5rocks.com/en/tutorials/speed/animations/
  latestKnownScrollY = window.scrollY;
  latestKnownY = $(window).height();
  requestTick();
}

function requestTick() {
  if (!ticking) {
    requestAnimationFrame(updateVisible);
  }
  ticking = true;
}

function updateVisible() {
  ticking = false; // Reset, ready to capture next scroll/resize
  chosenTree = localStorage.getItem("chosenTree");
  var currentScrollY = latestKnownScrollY;
  var currentHeight = latestKnownY;
  $("svg.tree").each(function(index, item) {
    if (onScreen($(this), currentScrollY, currentHeight)) {
      if ($(this).data('tree') != chosenTree) {
        var oldTree = $(item);
        $.ajax({
           url: "`r getOption('ProjectName')`_files/figure-html/tree" + chosenTree
                + "-char" + $(oldTree).data("char") + ".svg",
          context: item,
          datatype: 'xml'
        }).done(function (newTree) {
          $newSvg = $('svg', newTree);
          $newSvg.data('tree', chosenTree);
          $newSvg.click(highlightTaxon);
          $('text', $newSvg).each(applyHighlighting);
          $(this).replaceWith($newSvg);
        })
      }
    }
  });
  return false;
}

function switchTree(input) {
  newValue = input.value * 1;
  if (!isNaN(newValue)) {
    newValue = Math.floor(newValue);
    newValue = Math.max(newValue, 1);
    newValue = Math.min(newValue, input.max);
    $("input.switcherNumber").val(newValue);
    localStorage.setItem("chosenTree", newValue);
    updateVisible()
  }
  return false;
}


$(document).ready(function() {
  window.$my = {
    trees : $("svg.tree")
  };
  // Produce detail toggler for each character
  var originalToggler = $("div.toggleDetails");
  $stateNotes = $('div.state-notes');
  $stateNotes.each(function(index) {
      var $theseNotes = $(this);
      $theseNotes.attr('id', 'stateNotes' + index);
      var togglerClone = originalToggler
            .clone()
            .data('index', index)
            .click(toggleDetails)
            .insertBefore($theseNotes);
    }
  );
  $stateNotes.hide();
  originalToggler.remove();

  // Produce switcher box by each SVG plot
  var originalSwitcher = $("div.switcher");
  $my.trees.each(function() {
      var thisTree = $(this);
      var switcherClone = originalSwitcher.clone().insertBefore(thisTree);
      $('input.switcherNumber', switcherClone).data("char", thisTree.data("char"));
    }
  );

  localStorage.setItem("chosenTree", originalSwitcher.val());
  $my.trees.click(highlightTaxon)
  $my.trees.data('tree', originalSwitcher.val());
  $('text', $my.trees).each(applyHighlighting)

  $(".book-body").scroll(onScroll);
  $(".book-body").resize(onScroll);

  $(".body-inner").resize(onScroll);
  $(".body-inner").scroll(onScroll);
  // page-wrapper and page-inner are the other things that might trigger events?

  originalSwitcher.remove();
  return false;
});
