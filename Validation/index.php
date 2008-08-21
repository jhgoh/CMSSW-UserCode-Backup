<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN"
	"http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd">
<html xmlns="http://www.w3.org/1999/xhtml" xml:lang="en" lang="en">

<head>

<meta http-equiv="Content-Type" content="text/html;charset=utf-8" />
<title>Validation Plots</title>	 
<link rel="stylesheet" type="text/css" href="screen.css" />
<style type="text/css">@import url('screen.css');</style>

<?
	@include('/users/jhgoh/public_html/CMS/Validation/functions.php');

	#Initialize variables

	$CGIPATH = '/~jhgoh/python/rootView.py';
	$DATAPATH = '/users/jhgoh/public_html/CMS/Validation/data';

	$notify_msg = '';

	# Fill CMSSW_RELEASES by reading directory content under $DATADIR
	$releases_src = Array();
	if ( $dh = @opendir($DATAPATH) ){
		while ( False != ($f = @readdir($dh)) ){
			if ( $f == '.' or $f == '..' ) continue;
			if ( is_dir("$DATAPATH/$f") ) $releases_src[] = $f;
		}
		@closedir($dh);
	}
	else {
		$notify_msg .= "Alert : Cannot find data directory\n";
	}

	$validators_src = getRootObjects('TDirectoryFile', '/DQMData/([a-zA-Z0-9_]*V)$');

	$releases = Array();
	$subsystems = Array(); 

?>

<script type="text/javascript" language="javascript">
<!--
	function toggleHide(objId) {
		var obj = document.getElementById(objId);
		if ( obj.style.display == 'none' ) {
			obj.style.display = 'block';
		}
		else {
			obj.style.display = 'none';
		}
	}
-->
</script>

</head>

<body onload="">

<?	if ( $notify_msg != "" ) echo "<div id=\"notify\"><?=$notify_msg?></div>" ?>

<div id="header">
 <h1><a href="index.php">Validation Plots for CMSSW</a></h1>
</div>

<div id="main">
<form id="form" name="form" action="#" method="post">

 <div id="menu">
  <h2>Presets</h2>

  <h2>Data selection</h2>

  <h3><a href="#" onclick="toggleHide('div_validator')">Validator</a></h3>
  <div id="div_validator">
<?
	foreach ( $validators_src as $validator_src ) {
?>
   <input type="radio" name="validator" value="<?=$validator_src?>" onchange="this.form.submit()"<?=radio_postfix('validator', $validator_src)?>><?=$validator_src?></input><br/>
<?
	}

	if ( array_key_exists('validator', $_POST) ) {
		$validator = $_POST['validator'];
?>
  </div>

  <h3><a href="#" onclick="toggleHide('div_subsystem')">Subsystem</a></h3>
  <div id="div_subsystem">
<?
		$subsystems_src = getRootObjects('TDirectoryFile', "/DQMData/$validator/([A-Za-z0-9_]*)$");
		foreach ( $subsystems_src as $subsystem_src ) {
			if ( array_key_exists('subsystems', $_POST) and 
				 array_key_exists($subsystem_src, $_POST['subsystems']) ) {
				$subsystems[$subsystem_src] = Array();
				$postfix = ' checked="on"';
			}
			else {
				$postfix = '';
			}
			echo "<input type=\"checkbox\" name=\"subsystems[$subsystem_src]\"$postfix>$subsystem_src</input><br/>\n";

			$selectors_src = getRootObjects('TDirectoryFile', "DQMData/$validator/$subsystem_src/([A-Za-z0-9_]*)$");
			foreach ( $selectors_src as $selector_src ) {
				if ( array_key_exists('selector', $_POST) and
					 array_key_exists($subsystem_src, $_POST['selector']) and
					 array_key_exists($selector_src, $_POST['selector'][$subsystem_src]) ) {
					$subsystems[$subsystem_src][] = $selector_src;
					$postfix = ' checked="on"';
				}
				else {
					$postfix = '';
				}
				$short_name = substr($selector_src, 0, 20);
				echo "&nbsp;&nbsp;&nbsp;<input type=\"checkbox\" name=\"selector[$subsystem_src][$selector_src]\"$postfix>$short_name...</input><br/>\n";
			}
		}
	}
?>
  </div>

  <h3><a href="#" onclick="toggleHide('div_releases')">Releases</a></h3>
  <div id="div_releases">
<?
	foreach ( $releases_src as $release_src ) {
		if ( array_key_exists('releases', $_POST) and
			 array_key_exists($release_src, $_POST['releases']) ) {
			$releases[] = $release_src;
			$postfix = ' checked="checked"';
		}
		else { 
			$postfix = '';
		}
		echo "<input type=\"checkbox\" name=\"releases[$release_src]\"$postfix>$release_src</input><br/>\n";
	}
?>
  </div>

  <h3><a href="#" onclick="toggleHide('div_samples')">Samples</a></h3>
  <div id="div_samples">
<?
	$samples_src = Array('RelValSingleMuPt10', 'RelValSingleMuPt100', 'RelValSingleMuPt1000', 'RelValWM');
	foreach ( $samples_src as $sample_src ) {
		if ( array_key_exists('samples', $_POST) and 
			 array_key_exists($sample_src, $_POST['samples']) ) {
			$postfix = ' checked="on"';
		}
		else {
			$postfix = '';
		}
		echo "<input type=\"checkbox\" name=\"samples[$sample_src]\"$postfix>$sample_src</input><br/>\n";
	}
?>
  </div>

  <br/>

  <h2>Plot styles</h2>

  <label for="display_mode">Display mode</label>
  <select id="display_mode" name="display_mode">
   <option id="thumbnail" value="thumbnail"<?=select_postfix('display_mode', 'thumbnail')?>>Thumbnail</option>
   <option id="2column" value="2column"<?=select_postfix('display_mode', '2column')?>>Double column</option>
   <option id="1column" value="1column"<?=select_postfix('display_mode', '1column')?>>Single column</option>
  </select>
  <br/>

  <label for="collate_mode">Collate by : </label>
  <select id="collate_mode" name="collate_mode">
   <option value="none"<?=select_postfix('collate_mode', 'none')?>>None</option>
   <option value="release"<?=select_postfix('collate_mode', 'release')?>>Release</option>
   <option value="sample"<?=select_postfix('collate_mode', 'sample')?>>Sample</option>
  </select>
  <br/>

<!--
  <input type="checkbox" id="hide_stat" name="hide_stat"<?=checkbox_postfix('hide_stat')?>>Hide statistics</input><br/>

  <label for="palette">Color palette</label>
  <select id="palette" name="palette">
   <option id="Rainbow"<?=select_postfix('palette', 'Rainbow')?>>Rainbow</option>
  </select>
  <br/>
-->

  <label for="2d_style">2-D histogram style</label>
  <select id="2d_style" name="2d_style">
   <option id="COLZ" value="COLZ"<?=select_postfix('2d_style', 'COLZ')?>>COLZ</option>
   <option id="COL" value="COL"<?=select_postfix('2d_style', 'COL')?>>COL</option>
   <option id="LEGO" value="LEGO"<?=select_postfix('2d_style', 'LEGO')?>>LEGO</option>
  </select>
  <br/>

  <input type="submit" action="#" value="Apply"></input>
 </div>

 <div id="contents">
<?
	$validator = array_key_exists('validator', $_POST) ? $_POST['validator'] : '';

//	$subsystems_map = Array();
//	$selectors_map = Array();
//	if ( array_key_exists('subsystems', $_POST) ) {
//		$subsystems_map = $_POST['subsystems'];
//	}
	$samples_map = array_key_exists('samples', $_POST) ? $_POST['samples'] : Array();

	$display_mode = array_key_exists('display_mode', $_POST) ? $_POST['display_mode'] : '';
	$collate_mode = array_key_exists('collate_mode', $_POST) ? $_POST['collate_mode'] : '';

	$glb_opt_str = '';

	if ( $validator && $subsystems ) {
		foreach ( $subsystems as $subsystem=>$selectors ) {
			$histograms = getRootObjects('TH', "(/DQMData/$validator/$subsystem/[A-Za-z0-9_]*)$");
			if ( $collate_mode == 'release' ) {
				
				foreach ( $samples_map as $sample=>$onOff ) {
					echo "  <h3>$subsystem - $sample</h3>";
					foreach ( $histograms as $histogram ) {
						$histo_dir_paths = explode('/', $histogram);
						$last_idx = count($histo_dir_paths)-1;
						$title = $histo_dir_paths[$last_idx];

						$full_path_list = "";
						foreach ( $releases as $release ) {
							$full_path_list .= "$release/$sample.root:$histogram,";
						}
						$img_url = "$CGIPATH?h=$full_path_list";

						if ( $display_mode == 'thumbnail' ) showThumbnail($title, $img_url);
						else showDetail($title, $img_url, '-');

					}
					foreach ( $selectors as $selector ) {
						echo "   <h4>$selector</h4>";
						$selector_histograms = getRootObjects('TH', "(/DQMData/$validator/$subsystem/$selector/[A-Za-z0-9_]*)$");
						foreach ( $selector_histograms as $sel_histo ) {
							$histo_dir_paths = explode('/', $sel_histo);
							$last_idx = count($histo_dir_paths)-1;
							$title = $histo_dir_paths[$last_idx];

							$full_path_list = "";
							foreach ( $releases as $release ) {
								$full_path_list .= "$release/$sample.root:$sel_histo,";
							}
							$img_url = "$CGIPATH?h=$full_path_list";

							if ( $display_mode == 'thumbnail' ) showThumbnail($title, $img_url);
							else showDetail($title, $img_url, '-');
						}
					}
				}
			}
			else if ( $collate_mode == 'sample' ) {
				foreach ( $releases as $release ) {
					echo "  <h3>$subsystem - $release</h3>";
					foreach ( $histograms as $histogram ) {
						$histo_dir_paths = explode('/', $histogram);
						$last_idx = count($histo_dir_paths)-1;
						$title = $histo_dir_paths[$last_idx];

						$full_path_list = "";
						foreach ( $samples_map as $sample=>$onOff ) {
							$full_path_list .= "$release/$sample.root:$histogram,";
						}
						$img_url = "$CGIPATH?h=$full_path_list";

						if ( $display_mode == 'thumbnail' ) showThumbnail($title, $img_url);
						else showDetail($title, $img_url, '-');
					}
					foreach ( $selectors as $selector ) {
						echo "   <h4>$selector</h4>";
						$selector_histograms = getRootObjects('TH', "(/DQMData/$validator/$subsystem/$selector/[A-Za-z0-9_]*)$");
						foreach ( $selector_histograms as $sel_histo ) {
							$histo_dir_paths = explode('/', $sel_histo);
							$last_idx = count($histo_dir_paths)-1;
							$title = $histo_dir_paths[$last_idx];

							$full_path_list = "";
							foreach ( $samples_map as $sample=>$onOff ) {
								$full_path_list .= "$release/$sample.root:$sel_histo,";
							}
							$img_url = "$CGIPATH?h=$full_path_list";

							if ( $display_mode == 'thumbnail' ) showThumbnail($title, $img_url);
							else showDetail($title, $img_url, '-');
						}
					}
				}
			}
			else {
				foreach ( $releases as $release ) {
					foreach ( $samples_map as $sample=>$onOff ) {
						echo "  <h3>$subsystem - $release - $sample</h3>";
						foreach ( $histograms as $histogram ) {
							$histo_dir_paths = explode('/', $histogram);
							$last_idx = count($histo_dir_paths)-1;
							$title = $histo_dir_paths[$last_idx];

							$full_path_list = "$release/$sample.root:$histogram";
							$img_url = "$CGIPATH?h=$full_path_list";

							if ( $display_mode == 'thumbnail' ) showThumbnail($title, $img_url);
							else showDetail($title, $img_url, '-');
						}
						foreach ( $selectors as $selector ) {
							echo "   <h4>$selector</h4>";
							$selector_histograms = getRootObjects('TH', "(/DQMData/$validator/$subsystem/$selector/[A-Za-z0-9_]*)$");
							foreach ( $selector_histograms as $sel_histo ) {
								$histo_dir_paths = explode('/', $sel_histo);
								$last_idx = count($histo_dir_paths)-1;
								$title = $histo_dir_paths[$last_idx];

								$full_path_list = "$release/$sample.root:$sel_histo,";
								$img_url = "$CGIPATH?h=$full_path_list";

								if ( $display_mode == 'thumbnail' ) showThumbnail($title, $img_url);
								else showDetail($title, $img_url, '-');
							}
						}
					}
				}
			}
		}
	}
?>
 </div>

</form>
</div>

<div id="tail">
 <a href="http://validator.w3.org/check?uri=referer" target="_blank">XHTML Transitional 1.0</a> /
 <a href="http://jigsaw.w3.org/css-validator/check?uri=referer" target="_blank">CSS</a> / 
 Tested at Mozilla-Firefox<br/>
 Last updated : <?=date('Y M d', filemtime("index.php"))?><br/>
 Junghwan Goh (jhgoh@fnal.gov)
</div>

</body>

</html>
