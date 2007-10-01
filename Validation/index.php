<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN"
    "http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd">
<html xmlns="http://www.w3.org/1999/xhtml" xml:lang="en" lang="en">

<head>

<meta http-equiv="Content-Type" content="text/html;charset=utf-8" />
<title>Validation Plots</title>
<link rel="stylesheet" type="text/css" href="screen.css" />
<style type="text/css">@import url('screen.css');</style>

<script type="text/javascript">
</script>

</head>

<body>
<form id="subForm" action="<?=$PHP_SELF?>">
<?
	$webpath = '/afs/cern.ch/cms/Physics/muon/CMSSW/Performance/RecoMuon/Validation';

// Categorization of plots
//   add your validator and plots here.
	switch ($validator) {
	case 'RecoMuonValidator' :
		$categories = array('efficiency'=>array('GlbSim_effEta', 'StaSim_effEta','SeedSim_effEta','GlbSta_effEta','GlbSeed_effEta','StaSeed_effEta','GlbTk_effEta'),
				    'resolution'=>array('GlbEtaVsErrQPt_2', 'StaEtaVsErrQPt_2'));
		break;
	case 'MultiTrackAnalyzer':
		$categories = array('efficiency'=>array('effic','efficPt','fakes','fakerate'), 
				    'Pt resolution'=>array('ptres_vs_eta','ptres_vs_eta_1','ptres_vs_eta_2','ptres_vs_eta_chi2'),
				    'Pt pull'=>array('ptpull_vs_eta','ptpull_vs_eta_1','ptpull_vs_eta_2','ptpull_vs_eta_chi2'));
		break;
	}

	$cmssw_versions = array();
	$samples = array();
	$validators = array();
	$selectors = array();

	if ( $view != "twoColumn" ) $view = "thumbnail";
?>

<? 
printHead(); 
?>

<?
	if ( ($cmssw_versions = getListOfDir($webpath.'/data')) == array() ) {
		error("Directory $webpath is not ready");
		printTail();
		die();
	}
	
	if ( $cmssw_version != "" and is_dir("$webpath/data/$cmssw_version") ) {
		$samples = getListOfDir("$webpath/data/$cmssw_version");
		if ( $sample != "" and is_dir("$webpath/data/$cmssw_version/$sample") ) {
			$validators = getListOfDir("$webpath/data/$cmssw_version/$sample");
			if ( $validator != "" and is_dir("$webpath/data/$cmssw_version/$sample/$validator") ) {
				$selectors = getListOfDir("$webpath/data/$cmssw_version/$sample/$validator");
			}
		}
	}
?>
<div id="menu">
 <p>
  <label for="cmssw_version">Release :</label>
  <select name="cmssw_version" onchange="this.form.submit();">
   <option value="">=== CMSSW_VERSION ===</option> 
<?
	foreach ($cmssw_versions as $item) { 
?>
   <option id="<?=$item?>" value="<?=$item?>"<? if ($cmssw_version==$item) print ' selected="selected"'?>><?=$item?></option> 
<?
 	}  
?>
  </select>

<?
	if ( $cmssw_version != "" ) { 
?>
  <label for="sample">Sample :</label>
  <select name="sample" onchange="this.form.submit();">
   <option value="">=== Sample ===</option> 
<?
		foreach ($samples as $item) { 
?> 
   <option id="<?=$item?>" value="<?=$item?>"<? if ($sample==$item) print ' selected="selected"'?>><?=$item?></option> 
<?
		} 
?> 
  </select>
<?
	}

	if ( $cmssw_version != "" and $sample != "" ) {
?>
  <label for="validator">Validator :</label>
  <select name="validator" onchange="this.form.submit();">
   <option value="">=== Validator ===</option>
<?
		foreach ($validators as $item) { 
?>
   <option id="<?=$item?>" value="<?=$item?>"<? if ($validator==$item) print ' selected="selected"'?>><?=$item?></option>
<?
		} 
?>
  </select>
<?
	} 

	if ( $cmssw_version != "" and $sample != "" and $validator != "" ) {
?>
  <label for="selector">Selector :</label>
  <select name="selector" onchange="this.form.submit();">
   <option value="">=== Selector ===</option>
<?	    
		foreach ($selectors as $item)  {
?>
   <option id="<?=$item?>" value="<?=$item?>"<? if ($selector==$item) print ' selected="selected"'?>><?=$item?></option>
<?
		}
?>
</select>
    <? }

	if ( $cmssw_version != "" and $sample != "" and $validator != "" and $selector != "" ) {
?>
  <label for="category">Category :</label>
  <select name="category" onchange="this.form.submit();">
   <option value="">=== All Plots ===</option>
<?
		reset($categories);
		foreach (array_keys($categories) as $item)  {
?>
   <option id="<?=$item?>" value="<?=$item?>"<? if ($category==$item) print ' selected="selected"'?>><?=$item?></option>
<?
		}
?>
  </select>
<?
	}

?>
  <label for="keywords">Keyword :</label> <input type="text" name="keywords" value="<?=$keywords?>" class="text"/><br/>
  <label for="view">View</label><br/> 
   <input type="radio" name="view" value="thumbnail" <? if ($view == 'thumbnail') { ?> checked="checked"<? } ?> /> thumbnail<br/>
   <input type="radio" name="view" value="twoColumn" <? if ($view == 'twoColumn') { ?> checked="checked"<? } ?> /> twoColumn<br/>
<!--<input type="checkbox" name="thumbnail" <? if ($thumbnail != "") { ?>checked="checked"<? } ?> /> <br/>-->
  <input type="hidden" name="mode" value="<?=$mode?>"/>
  <input type="submit" name="OK" value="OK"/>
 </p>
</div>

<div id="content">
<?
	// Print validation plots if we set CMSSW versions correctly
	if ( $sample == "" or $validator == "" ) {
		if ( $cmssw_version != "" ) {
?>
 <h2>Release <span style="color:red"><?=$cmssw_version?></span></h2>
<?
			$desc = "$webpath/data/$cmssw_version/description.txt";
			if ( is_file($desc) ) readfile($desc);
			print "</div>";
			printTail();
			die();
		}
		else {
			error("Please choose menu items");
			print "</div>";
			printTail();
			die();
		}
	}
?>
 <h2>Validation plots of <span style="color:red"><?=$cmssw_version?></span>, <span style="color:red"><?=$sample?></span></h2>
 <h3>Summary</h3>
 <pre class="description">
<?
	$workArea = "$webpath/data/$cmssw_version/$analyzer/$sample/$validator/$selector";
	if ( is_file("$workArea/description.txt") ) readfile("$workArea/description.txt");
?>
 </pre>

<? 
	if ( $category == "" ) $images = getListOfImages($workArea);
	else $images = $categories[$category];

	$images = filter($images, explode(' ', $keywords));
?>
 <h3>Result (<?=count($images)?> plots)</h3>
<?
	if ( $view == "thumbnail" ) {
		foreach ($images as $item) {
			$img  = "data/$cmssw_version/$sample/$validator/$selector/$item.gif";
			$link = "index.php?cmssw_version=$cmssw_version&amp;sample=$sample&amp;validator=$validator&amp;selector=$selector&amp;keywords=$item&amp;view=twoColumn";
?>
 <table class="thumbnail">
  <tr><th><?=$item?></th></tr>
  <tr><td>
    <a href="<?=$link?>"><img src="<?=$img?>" alt="<?=$item?>" /></a>
  </td></tr>
 </table>
<?
		}
	}
	else {
?>
 <table class="fullTable">
<? 
		if ( count($images) == 1 ) {
			$item = $images[0];
			$img  = "data/$cmssw_version/$sample/$validator/$selector/$item.gif";
?>
  <tr><th><?=$item?></th></tr>
  <tr><td><img src="<?=$img?>" alt="<?=$item?>"/></td></tr>
  <tr><td><p><?if (is_file("$workArea/$item.txt")) readfile("$workArea/$item.txt")?></p></td></tr>
<?
		}
		elseif ( count($images) >= 2 ) {
			$n = 2*floor(count($images)/2);
			for ($i=0; $i<$n; $i+=2) {
				$item1 = $images[$i];
				$item2 = $images[$i+1];
				$img1  = "data/$cmssw_version/$sample/$validator/$selector/$item1.gif";
				$img2  = "data/$cmssw_version/$sample/$validator/$selector/$item2.gif";
				$link1 = "index.php?cmssw_version=$cmssw_version&amp;sample=$sample&amp;validator=$validator&amp;selector=$selector&amp;keywords=$item1";
				$link2 = "index.php?cmssw_version=$cmssw_version&amp;sample=$sample&amp;validator=$validator&amp;selector=$selector&amp;keywords=$item2";
?>
   <tr><th><?=$item1?></th><th><?=$item2?></th></tr>
   <tr><td><a href="<?=$link1?>"><img src="<?=$img1?>" alt="<?=$item1?>"/></a></td>
       <td><a href="<?=$link2?>"><img src="<?=$img2?>" alt="<?=$item2?>"/></a></td></tr>
   <tr><td><pre><?if (is_file("$workArea/$item1.txt")) readfile("$workArea/$item1.txt")?></pre></td>
       <td><pre><?if (is_file("$workArea/$item2.txt")) readfile("$workArea/$item2.txt")?></pre></td></tr>
<? 
			}
			if ($n != count($images) ) {
				$item1 = $images[$n];
				$img1  = "data/$cmssw_version/$sample/$validator/$selector/$item1.gif";
				$link1 = "index.php?cmssw_version=$cmssw_version&amp;sample=$sample&amp;validator=$validator&amp;selector=$selector&amp;keywords=$item1";
?>
   <tr><th><?=$item1?></th><th>&nbsp;</th></tr>
   <tr><td><a href="<?=$link1?>"><img src="<?=$img1?>" alt="<?=$item1?>"/></a></td>
       <td>&nbsp;</td></tr>
   <tr><td><pre><?if (is_file("$workArea/$item1.txt")) readfile("$workArea/$item1.txt")?></pre></td>
       <td>&nbsp;</td></tr>
<?
			}
		}
?>
 </table>
<? 
	}
?>
</div>
<?
	printTail(); 
?>

<?//////////////////////////////////Functions///////////////////////?>

<?
	function error($msg) {
?>
<div id="error">
 <h2>Error!!!</h2>
 <p><?=$msg?></p>
</div>
<?
	} 

	function getListOfDir($path) {
		$dirList = array();
		if ( $dh = @opendir($path) ) {
			while ( False !== ($f = @readdir($dh) ) ) {
				if ( $f != '.' and $f != '..' ) array_push($dirList, $f);
			}
			@closedir($dh);
		}
		return $dirList;
	}

	function getListOfImages($path) {
		$gifList = array();
		if ( $dh = @opendir($path) ) {
			while ( False !== ($f = @readdir($dh) ) ) {
				if ( strlen($f) < 5 ) continue;
				if ( substr($f, -4) != ".gif") continue;
				array_push($gifList, substr($f, 0, -4));
			}
			@closedir($dh);
		}
		return $gifList;
	}

	function filter($strings, $keywords) {
		$strList = array();
		if ( count($strings) == 0 or count($keywords) == 0 ) return $strings;
		foreach ($strings as $str) {
			$matchedAll = True;
			foreach ($keywords as $keyword) {
				if ( $keyword == "" ) continue;
				if ( ! stristr($str, $keyword) ) $matchedAll = False;
			}
			if ( $matchedAll ) array_push($strList, $str);
		}
		return $strList;
	}

	function printHead() { ?>
<div id="header">
 <!--<img src="http://cmsdbs.cern.ch/images/CMSLogo.gif" style="width:30px;height:30px;float:left;padding-right:15px;" alt="CMS logo"/>-->
 <h1><a href="index.php">Validation Plots for CMSSW</a></h1>
</div>

<div id="main">
<? 
	} 

	function printTail() { ?>
</div>

<div id="tail">
 <a href="http://validator.w3.org/check?uri=referer" target="_blank">XHTML Transitional 1.0</a> /
 <a href="http://jigsaw.w3.org/css-validator/check?uri=referer" target="_blank">CSS</a> / 
 Tested at Mozilla-Firefox<br/>
 Last updated : 2007.Oct.01<br/>
 Junghwan Goh (jhgoh@fnal.gov)
</div>

</form>
</body>

</html>
<? 
	} 
?>

