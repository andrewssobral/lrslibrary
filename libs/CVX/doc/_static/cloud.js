/* ~~~~~~~~~~~~~~
 * cloud.js_t
 * ~~~~~~~~~~~~~~
 *
 * Various bits of javascript driving the moving parts behind various
 * parts of the cloud theme. Handles things such as toggleable sections,
 * collapsing the sidebar, etc.
 *
 * :copyright: Copyright 2011 by Assurance Technologies
 *             (Modified by Michael C. Grant, CVX Research, Inc.)
 * :license: BSD
 */

/*
 * MathJax fix
 * This reduces the height of each MathJax display element to a multiple of 20px for the purposes of vertical rhythm preservation. I should honestly change the hard-coded 20 to retrieve the local line_height. Oh well.
 */

MathJax.Hub.Register.MessageHook("End Process",function(message){
        $("div.MathJax_Display").each(function(){
                var trueHeight = $(this).height();
                var remainder = trueHeight%20;
                $(this).css({height:(trueHeight-remainder)+'px',position:'relative',top:-(remainder/2)+'px'});
       });
});

/* ==========================================================================
 * highlighter #2
 * ==========================================================================
 *
 * Sphinx's highlighter marks some objects when user follows link,
 * but doesn't include section names, etc. This catches those.
 */
$(document).ready(function (){
  // helper to locate highlight target based on #fragment
  function locate_target(){
    // find id referenced by #fragment
    var hash = document.location.hash;
    if(!hash) return null;
    var section = document.getElementById(hash.substr(1));
    if(!section) return null;

    // could be div.section, or hidden span at top of div.section
    var name = section.nodeName.toLowerCase();
    if(name != "div"){
      if(name == "span" && section.innerHTML == "" &&
         section.parentNode.nodeName.toLowerCase() == "div"){
        section = section.parentNode;
      }
    }
    // now at section div and either way we have to find title element - h2, h3, etc.
    var children = $(section).children("h2, h3, h4, h5, h6");
    return children.length ? children : null;
  }

  // init highlight
  var target = locate_target();
  if(target) target.addClass("highlighted");

  // update highlight if hash changes
  $(window).bind("hashchange", function () {
    if(target) target.removeClass("highlighted");
    target = locate_target();
    if(target) target.addClass("highlighted");
  });
});

/* ==========================================================================
 * toggleable sections
 * ==========================================================================
 *
 * Added expand/collapse button to any collapsible RST sections.
 * Looks for sections with CSS class "html-toggle",
 * along with the optional classes "expanded" or "collapsed".
 * Button toggles "html-toggle.expanded/collapsed" classes,
 * and relies on CSS to do the rest of the job displaying them as appropriate.
 */

$(document).ready(function (){
  function init(){
    // get header & section, and add static classes
    var header = $(this);
    var section = header.parent();
    header.addClass("html-toggle-button");

    // helper to test if url hash is within this section
    function contains_hash(){
      var hash = document.location.hash;
      return hash && (section[0].id == hash.substr(1) ||
              section.find(hash.replace(".","\\.")).length>0);
    }

    // helper to control toggle state
    function set_state(expanded){
      if(expanded){
        section.addClass("expanded").removeClass("collapsed");
        section.children().show();
      }else{
        section.addClass("collapsed").removeClass("expanded");
        section.children().hide();
        section.children("span:first-child:empty").show(); /* for :ref: span tag */
        header.show();
      }
    }

    // initialize state
    set_state(section.hasClass("expanded") || contains_hash());

    // bind toggle callback
    header.click(function (){
      set_state(!section.hasClass("expanded"));
      $(window).trigger('cloud-section-toggled', section[0]);
    });

    // open section if user jumps to it from w/in page
    $(window).bind("hashchange", function () {
      if(contains_hash()) set_state(true);
    });
  }

  $(".html-toggle.section > h2, .html-toggle.section > h3, .html-toggle.section > h4, .html-toggle.section > h5, .html-toggle.section > h6").each(init);
});
/* ==========================================================================
 * collapsible sidebar
 * ==========================================================================
 *
 * Adds button for collapsing & expanding sidebar,
 * which toggles "document.collapsed-sidebar" CSS class,
 * and relies on CSS for actual styling of visible & hidden sidebars.
 */

$(document).ready(function (){
  var holder = $('<div class="sidebartoggle"><button id="sidebar-hide" title="click to hide the sidebar">&laquo;</button><button id="sidebar-show" style="display: none" title="click to show the sidebar">sidebar &raquo;</button></div>');
  var doc = $('div.document');

  var show_btn = $('#sidebar-show', holder);
  var hide_btn = $('#sidebar-hide', holder);
  var copts = { expires: 7, path: DOCUMENTATION_OPTIONS.url_root };

  show_btn.click(function (){
    doc.removeClass("collapsed-sidebar");
    hide_btn.show();
    show_btn.hide();
    $.cookie("sidebar", "expanded", copts);
    $(window).trigger("cloud-sidebar-toggled", false);
  });

  hide_btn.click(function (){
    doc.addClass("collapsed-sidebar");
    show_btn.show();
    hide_btn.hide();
    $.cookie("sidebar", "collapsed", copts);
    $(window).trigger("cloud-sidebar-toggled", true);
  });

  var state = $.cookie("sidebar");


  doc.append(holder);

  if (state == "collapsed"){
    doc.addClass("collapsed-sidebar");
    show_btn.show();
    hide_btn.hide();
  }
});

/* ==========================================================================
 * sidebar toc highlighter
 * ==========================================================================
 *
 * highlights toc entry for current section being viewed.
 * only enabled under stick mode
 */
$(document).ready(function (){
  // pre-lookup all the links, hack in a class for css styling.
  var links = $("div.sphinxsidebar p.logo + h3 + ul").addClass("sphinxtoclist").find("a");
  var h1_section = $("h1").parent();
  if ( $("div.toctree-wrapper").length == 0 ) {
      $("div.sphinxsidebar p.logo").remove();
  }

  // function to update toc markers
  function update_current(){
    // determine viewable range
    var start = $(window).scrollTop();
    var stop = start + $(window).height();

    // clear flags from last pass
    links.removeClass("toggled").removeClass("current");
    var link_count = links.length;

    // set 'current' class for all currently viewable sections in toc
    for(var i=0; i < link_count; ++i){
      var elem = $(links[i]); // XXX: could cache this in another list

      // hack to skip elements hidden w/in a toggled section
      if(elem.hasClass("toggled")) continue;

      // lookup section referenced by link
      var tag = elem.attr("href");
      var section = (tag == "#") ? h1_section : $(tag);

      // hack to flag subsections w/in a toggled section
      var toggle_children = null;
      if(section.is(".html-toggle.collapsed")){
        toggle_children = elem.parent().find("ul a").addClass("toggled");
      }

      // if section is off-screen, don't mark it
      var top = section.offset().top;
      if(top > stop || top + section.height() < start) continue;

      // if section has children, only count it if prologue is visible.
      var child = section.find("div.section"); // FIXME: only need first match
      if(child.length && child.offset().top < start) continue;

      // mark it!
      elem.addClass("current");
//      if(toggle_children) toggle_children.addClass("current");
    }
  }

  // run function now, and every time window is resized
  update_current();
  $(window).scroll(update_current)
           .resize(update_current)
           .bind('cloud-section-toggled', update_current)
           .bind('cloud-sidebar-toggled', update_current);
});


/* ==========================================================================
 * header breaker
 * ==========================================================================
 *
 * attempts to intelligently insert linebreaks into page titles, where possible.
 * currently only handles titles such as "module - description",
 * adding a break after the "-".
 */
$(document).ready(function (){
  // get header's content, insert linebreaks
  var header = $("h1");
  var orig = header[0].innerHTML;
  var shorter = orig;
  if($("h1 > a:first > tt > span.pre").length > 0){
      shorter = orig.replace(/(<\/tt><\/a>\s*[-:]\s+)/im, "$1<br> ");
  }
  else if($("h1 > tt.literal:first").length > 0){
      shorter = orig.replace(/(<\/tt>\s*[-:]\s+)/im, "$1<br> ");
  }
  if(shorter == orig){
    return;
  }

  // hack to determine full width of header
  header.css({whiteSpace: "nowrap", position:"absolute"});
  var header_width = header.width();
  header.css({whiteSpace: "", position: ""});

  // func to insert linebreaks when needed
  function layout_header(){
    header[0].innerHTML = (header_width > header.parent().width()) ? shorter : orig;
  }

  // run function now, and every time window is resized
  layout_header();
  $(window).resize(layout_header)
           .bind('cloud-sidebar-toggled', layout_header);
});

