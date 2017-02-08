/* Heavily edited version of Matt Kruse's list tree */
/* Feel free to use without crediting me, but please credit Matt */
/* WWW: http://www.mattkruse.com/ */

var treeC     = "mktree";
var nActiveC  = "liActive";
var nClosedC  = "liActive liClosed";
var nOpenC    = "liActive liOpen";
var nBulletC  = "liBullet";
var nLinkC    = "bullet";
var jsOnlyC   = "jsonly";
var jsHideC   = "jshide";
var jsMenuC   = "jsmenu";
var jsIFrameC = "jsframe";
function expandTree( treeId ) {
    var arr = document.getElementById( treeId ).getElementsByClassName( nClosedC );
    while ( arr.length ) arr[0].className = nOpenC;
}
function collapseTree( treeId ) {
    var arr = document.getElementById( treeId ).getElementsByClassName( nOpenC );
    while ( arr.length ) arr[0].className = nClosedC;
}
var convertTrees = function() {
    if ( !document.createElement ) 
        return;
    if ( window.location == window.parent.location )
		document.getElementsByTagName('body')[0].className = "no_iframe";       
	function scanList( arr ) {
		var retval = false;
		for ( var i = 0 ; i < arr.length ; ++i ) {
			var it = arr[i];
			switch ( it.nodeName ) {
			case "UL":
				scanList( it.childNodes );
				retval = true;
				break;
			case "LI": {
	            var s = document.createElement("SPAN");
        	    s.className = nLinkC;
				if ( scanList( it.childNodes ) ) {
					it.className = nClosedC;
					var fc = it.firstChild;
					switch ( fc.nodeName ) {
						case "#text": case "P": case "B": case "EM":
							s.appendChild( fc );
							break;
						default:
							s.appendChild( document.createTextNode( '\u00A0' ) );
					}
					s.onclick = function () {
						this.parentNode.className = ( this.parentNode.className == nOpenC ) ? nClosedC : nOpenC;
						return false;
					}
                } else
					it.className = nBulletC;
    	        it.insertBefore( s, it.firstChild );
			} break; }
		}
		return retval;
	}
    var els = document.getElementsByClassName( jsMenuC );
    for ( var i = 0 ; i < els.length ; ++i ) els[i].style.display = "inline";
    var els = document.getElementsByClassName( jsOnlyC );
    for ( var i = 0 ; i < els.length ; ++i ) els[i].style.display = "block";
    var els = document.getElementsByClassName( jsHideC );
    for ( var i = 0 ; i < els.length ; ++i ) els[i].style.display = "none";
    // var els = document.getElementsByClassName( jsIFrameC );
    // var ival = ( window.location != window.parent.location ) ? "none" : "inline";
    // for ( var i = 0 ; i < els.length ; ++i ) els[i].style.display = ival;
    scanList( document.getElementsByClassName( treeC ) );
	var cookie = document.cookie;
	cookie = cookie.match( '(^|;)mktree=([+-]*)(;|$)' );
	if ( cookie ) {
		cookie = cookie[2];
		var arr = document.getElementsByClassName( nActiveC );
		if ( arr.length == cookie.length )
			for ( var i = 0 ; i < arr.length ; ++i )
				arr[i].className = cookie[i] == '-' ? nOpenC : nClosedC;
	}				
}
var saveTrees = function() {
	var cookie = "mktree=";
	var arr = document.getElementsByClassName( nActiveC );
	for ( var i = 0 ; i < arr.length ; ++i )
		cookie += arr[i].className == nOpenC ? "-" : "+";
	cookie += ";"		
	document.cookie = cookie;
}
var q = function(a,b) { return a?(b?function(){a();b();}:a):b; };
window.onload = q(window.onload,convertTrees);
window.onunload = q(window.onunload,saveTrees);


