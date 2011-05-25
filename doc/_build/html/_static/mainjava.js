
var popupHandle;

function popup(picUrlString, windowWidth, windowHeight)
{

  // always close the old one , so only one at a time is open
  if(popupHandle || popupHandle!=null)
  {
    if (!popupHandle.closed) popupHandle.close();
  }
  popupHandle=null;

  // create a feature string for the popup
  var x=(screen.width-windowWidth)/2
  var y=(screen.height-windowHeight)/2
  var featureString = "toolbar=no,scrollbars=no,resizable=no"
  featureString = ',left='+x + ',top='+y
  featureString += ',width='+windowWidth+',height='+windowHeight

  // open the popup
  // We use a php script to work arround some browser specific problems
  // Opera doesn't seem to handle document.close() correctly on the popup.
  popupHandle = window.open( "popup.php3?windowWidth="+windowWidth+"&windowHeight="+windowHeight+"&picUrlString="+picUrlString ,"popup",featureString)
  return popupHandle;

}

function winclose()
{
  if (window.popupHandle!=null && !window.popupHandle.closed)
  {
    window.popupHandle.close();
    }
}

function doNothing(){} // does nothing but required by JavaScript in this case


