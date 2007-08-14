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
<?
	$webpath = '/afs/cern.ch/user/j/jhgoh/public/html/Validation';
	
	// Categorization of plots
	//   add your validator and plots here.
	switch ($validator) {
	case 'SingleMuonValidator' :
		$categories = array('efficiency'=>array('GlbSim_effEta', 'StaSim_effEta'),
				    'resolution'=>array('GlbEtaVsErrQPt_2', 'StaEtaVsErrQPt_2'));
		break;
	case 'MultiTrackAnalyzer':
		$categories = array('efficiency'=>array('eff'), 
				    'resolution'=>array('QPtResMean'));
		break;
	}

	$cmssw_versions = array();
	$samples = array();
	$validators = array();
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
		}
	}
?>
<div id="menu">
 <p>
  Release :
  <select name="cmssw_version" style="width:150px;">
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
  Sample : 
  <select name="sample" style="width:150px">
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
  Validator :
  <select name="validator" style="width:150px">
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
  Category :
  <select name="category" style="width:150px">
   <option value="">=== Category ===</option>
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
  Keyword : <input type="text" name="keywords" value="<?=$keywords?>"/><br/>
  Thumbnail <input type="checkbox" name="thumbnail" <? if ($thumbnail != "") { ?>checked="checked"<? } ?> /> <br/>
  <input type="hidden" name="mode" value="<?=$mode?>"/>
  <input type="submit" name="submit" value="OK"/>
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
	$workArea = "$webpath/data/$cmssw_version/$analyzer/$sample/$validator";
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
	if ( $thumbnail != "" ) {
		foreach ($images as $item) {
			$img  = "data/$cmssw_version/$sample/$validator/$item.gif";
			$link = "index.php?cmssw_version=$cmssw_version&sample=$sample&validator=$validator&keywords=$item";
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
			$img  = "data/$cmssw_version/$sample/$validator/$item.gif";
?>
  <tr><th><?=$item?></th></tr>
  <tr><td><img src="<?=$img?>" alt="<?=$item?>"/></td></tr>
  <tr><td><p><?if (is_file("$workArea/$item.txt")) readfile("$workArea/$item.txt")?></p></td></tr>
<?
		}
		elseif ( count($images) >= 2 ) {
			$n = 2*((count($images)+1)/2);
			for ($i=0; $i<$n; $i+=2) {
				$item1 = $images[$i];
				$item2 = $images[$i+1];
				$img1  = "data/$cmssw_version/$sample/$validator/$item1.gif";
				$img2  = "data/$cmssw_version/$sample/$validator/$item2.gif";
				$link1 = "index.php?cmssw_version=$cmssw_version&sample=$sample&validator=$validator&keywords=$item1";
				$link2 = "index.php?cmssw_version=$cmssw_version&sample=$sample&validator=$validator&keywords=$item2";
?>
   <tr><th><?=$item1?></th><th><?=$item2?></th></tr>
   <tr><td><a href="<?=$link1?>"><img src="<?=$img1?>" alt="<?=$item1?>"/></a></td>
       <td><a href="<?=$link2?>"><img src="<?=$img2?>" alt="<?=$item2?>"/></a></td></tr>
   <tr><td><pre><?if (is_file("$workArea/$item1.txt")) readfile("$workArea/$item1.txt")?></pre></td>
       <td><pre><?if (is_file("$workArea/$item2.txt")) readfile("$workArea/$item2.txt")?></pre></td></tr>
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
<form action="index.php">
<? 
	} 

	function printTail() { ?>
</form>
</div>

<div id="tail">
 <a href="http://validator.w3.org/check?uri=referer" target="_blank">XHTML Transitional 1.0</a> /
 <a href="http://jigsaw.w3.org/css-validator/check?uri=referer" target="_blank">CSS</a> / 
 Tested at Mozilla-Firefox<br/>
 Last updated : 2007.Aug.9<br/>
 Junghwan Goh (jhgoh@fnal.gov)
</div>

</body>

</html>
<? 
	} 
?>

