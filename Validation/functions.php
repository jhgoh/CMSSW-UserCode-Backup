<?php
	#Read list of validators from the *.root.ls file#
	function getRootObjects($type_match, $path_match)
	{
		$objects = Array();

		global $DATAPATH;
		if($fh = fopen("$DATAPATH/Validation.root.ls", 'r')){
			while ( !feof($fh)){
				$li = explode("\t", trim(fgets($fh)), 2);
				if(count($li) == 2)list($type, $path) = $li;
				if(!ereg($type_match, $type))continue;
				if(ereg($path_match, $path, $object_match)){
					$objects[] = $object_match[1];
				}
			}
			@fclose($fh);
		}

		return $objects;
	}

	function select_postfix($var_name, $option_name)
	{
		if(array_key_exists($var_name, $_POST) and
			 $_POST[$var_name] == $option_name){
			return ' selected="selected"';
		}
	}

	function checkbox_postfix($var_name)
	{
		if(array_key_exists($var_name, $_POST) and
			 $_POST[$var_name] != ''){
			return ' checked="checked"';
		}
	}

	function radio_postfix($var_name, $value)
	{
		if(array_key_exists($var_name, $_POST) and
			$_POST[$var_name] == $value){
			return ' checked="checked"';
		}
	}

	function showThumbnail($title, $img_url) {
?>
  <table class="thumbnail">
   <tr><th><?="$title"?></th></tr>
   <tr><td><a href="<?="$img_url&height=500&width=600"?>" target="_blank"><img src="<?="$img_url&width=240&height=200"?>"/></a></td></tr>
  </table>
<?
	}

	function showDetail($title, $img_url, $desc) {
?>
  <table class="detail">
   <tr><th colspan="2"><?=$title?></th></tr>
   <tr><td><a href="<?="$img_url&height=500&width=600"?>" target="_blank"><img src="<?=$img_url?>"/></a></td>
       <td><?=$desc?></td></tr>
  </table>
<?
	}


?>
