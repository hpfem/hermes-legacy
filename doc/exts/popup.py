from docutils import nodes
from docutils.parsers.rst import directives

CODE = """
<script type="text/javascript">
<!--
var popupHandle;

function popup(myPopup)
{
  var win = "videos.html" + myPopup;
  // always close the old one , so only one at a time is open
  if(popupHandle || popupHandle!=null)
  {
    if (!popupHandle.closed) popupHandle.close();
  }
  popupHandle=null;
   if(true)
 {
  // open the popup
  popupHandle = window.open( win, "myWindow", 
"height = 550, width = 900, resizable = 0" )
  return popupHandle;
 }
}

function winclose()
{
  if (window.popupHandle!=null && !window.popupHandle.closed)
  {
    window.popupHandle.close();
    }
}

function doNothing(){} // does nothing but required by JavaScript in this case
//-->
</script>

<p><b><a href="javascript:doNothing()" onClick="popupHandle=popup(%(name)s)">Watch the Video</a></b></p>

"""

PARAM = """\n    <param name="%s" value="%s"></param>"""

def popup(name, args, options, content, lineno,
            contentOffset, blockText, state, stateMachine):
    """ Restructured text extension for popups """
    if len(content) == 0:
        return
    string_vars = {
        'name': content[0],
        'extra': ''
        }
    extra_args = content[1:] # Because content[0] is ID
    extra_args = [ea.strip().split("=") for ea in extra_args] # key=value
    extra_args = [ea for ea in extra_args if len(ea) == 2] # drop bad lines
    extra_args = dict(extra_args)
    if extra_args:
        params = [PARAM % (key, extra_args[key]) for key in extra_args]
        string_vars['extra'] = "".join(params)
    return [nodes.raw('', CODE % (string_vars), format='html')]
popup.content = True
directives.register_directive('popup', popup)
