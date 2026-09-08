// Historical v5 native chart layout. Inputs are supplied by the caller.
import fs from 'node:fs/promises';
export async function addHistoricalGenomicControl({presentation, simulation:d, equations:math}) {
const C={ink:'#16313F',gray:'#53656F',light:'#E5EBEE',blue:'#0072B2',orange:'#E69F00'};
const p=presentation;
const s=p.slides.add(); s.background.fill='#FFFFFF';
function text(t,x,y,w,h,size=24,opts={}) {
 const o=s.shapes.add({name:opts.name||t,geometry:'textbox',position:{left:x,top:y,width:w,height:h},fill:'none',line:{fill:'none',width:0}});
 o.text=t; o.text.style={typeface:'Calibri',fontSize:size,color:opts.color||C.ink,bold:opts.bold||false,alignment:opts.align||'left',verticalAlignment:'top',autoFit:'none',wrap:'square',insets:{left:0,right:0,top:0,bottom:0}};
 return o;
}
function line(x1,y1,x2,y2,color=C.gray,width=2) {
 const x=Math.min(x1,x2),y=Math.min(y1,y2),w=Math.max(Math.abs(x2-x1),.1),h=Math.max(Math.abs(y2-y1),.1);
 return s.shapes.add({name:'GC connector',geometry:'custom',position:{left:x,top:y,width:w,height:h},fill:'none',line:{fill:color,width},customPaths:[{width:w,height:h,commands:[{moveTo:{x:x1-x,y:y1-y}},{lineTo:{x:x2-x,y:y2-y}}]}]});
}
function arrow(x1,y1,x2,y2) {line(x1,y1,x2,y2,C.gray,2.5);line(x2-11,y2-7,x2,y2,C.gray,2.5);line(x2-11,y2+7,x2,y2,C.gray,2.5);}
async function equation(key,x,y,width,maxHeight=70) {
 const e=math[key], scale=Math.min(width/e.width,maxHeight/e.height);
 s.images.add({blob:new Uint8Array(await fs.readFile(e.file)),contentType:'image/svg+xml',fit:'contain',alt:(e.number?'Equation '+e.number+'. ':'Panel annotation. ')+e.latex,position:{left:x,top:y,width:e.width*scale,height:e.height*scale}});
}
const axis=(title)=>({title:{text:title,textStyle:{typeface:'Calibri',fontSize:21,fill:C.ink}},min:0,max:32,majorUnit:10,numberFormatCode:'0',textStyle:{typeface:'Calibri',fontSize:19,fill:C.gray},line:{fill:C.gray,width:1},majorGridlines:{fill:C.light,width:.7},minorGridlines:null});
function panel(title,x,yValues,color,yTitle) {
 s.charts.add('scatter',{
  position:{left:x,top:262,width:492,height:331},title,titlePlacement:'aboveChart',
  titleTextStyle:{typeface:'Calibri',fontSize:27,fill:color,bold:true},
  series:[
   {name:'Null expectation',xValues:[0,32],values:[0,32],fill:C.gray,line:{fill:C.gray,width:1.7},marker:{symbol:'none'}},
   {name:'Simulated null SNPs',xValues:d.coordinates.expected,values:yValues,fill:color,line:{fill:'none',width:0},marker:{symbol:'circle',size:3}}
  ],scatterOptions:{style:'lineWithMarkers'},hasLegend:false,
  xAxis:axis('Expected χ²₁ quantile'),yAxis:axis(yTitle),
  chartFill:'#FFFFFF',chartLine:{fill:'none',width:0},plotAreaFill:'#FFFFFF',plotAreaLine:{fill:'none',width:0}
 });
}
text('DIAGNOSING INFLATION',80,17,1120,23,15,{name:'Keep source section',color:C.gray});
text('Genomic control',80,43,1120,99,44,{name:'Keep source title',bold:true});
text('Devlin & Roeder (1999), Biometrics;  Bulik-Sullivan et al. (2015a), Nature Genetics',80,687,1120,25,18,{name:'Keep source citation',color:C.gray});
text('68',1210,688,50,23,16,{name:'Keep source page',color:C.gray});
line(80,145,1200,145,C.light,2);
text('100,000 simulated null SNPs · uniform inflation',80,164,1120,37,27);
text('Estimate inflation',80,217,295,36,26,{bold:true});
await equation('estimate',390,201,540,55);
text('(20)',1158,214,46,32,24,{align:'right'});
panel('Before correction',70,d.coordinates.before,C.orange,'Observed test statistic');
panel('After correction',718,d.coordinates.after,C.blue,'Corrected test statistic');
text('Divide every\nstatistic by\n1.391',563,327,152,88,23,{align:'center',bold:true});
await equation('rescale',565,428,154,65);
text('(21)',616,494,48,30,22,{align:'center'});
arrow(579,535,700,535);
await equation('before_label',231,601,190,29);
await equation('after_label',880,601,190,29);
line(522,615,552,615,C.gray,1.7);
text('Null expectation',562,601,194,33,21,{color:C.gray});
text('Figure 26. Genomic control rescales test statistics; fitted SNP effects are unchanged.',80,650,1120,31,22,{color:C.gray,name:'Figure caption'});
s.speakerNotes.textFrame.setText(`Figure 26 is an original computed teaching simulation, not a published figure. ${d.notes}\n\n${d.n} independent z_i ~ N(0,1); T_i = 1.4 z_i^2. R seed 904, after the same runif draws as the subsequent LDSC example. ${d.R}. Generating inflation 1.4; estimated lambda_GC = ${d.lambda_before}. All statistics are divided by that same estimate. Median before = ${d.median_before}; after = ${d.median_after}; null median = ${d.null_median}. The corrected lambda equals 1 by construction, not by an independent calibration test. False-positive rate at nominal .05: ${d.false_positive_rate_0_05_before} before and ${d.false_positive_rate_0_05_after} after. All 100,000 statistics determine the estimate; 540 quantile positions are displayed, including all 200 largest statistics. Both panels show the same ranked SNPs with identical axes. No fitted effect estimates are modified.\n\nNotation matches Equations 19–21 in the lecture: T_i is the one-df test statistic; lambda_GC is the genomic-control inflation factor; T_i^GC is the rescaled statistic.\n\nReferences describe the method and its limitations, not the source of this synthetic figure:\nDevlin & Roeder (1999), Biometrics 55, 997–1004. https://doi.org/10.1111/j.0006-341x.1999.00997.x\nBulik-Sullivan et al. (2015a), Nature Genetics 47, 291–295. https://doi.org/10.1038/ng.3211`);
return s;
}
