
var mainPanel = ui.Panel({
  style: {width: '500px'}
});
// var c = ee.Number(0); c=0;

var title = ui.Label({
  value: 'Water occurrence in IGNP Command',
  style: {'fontSize': '20px',color:'blue',margin: '8px 0 -3px 8px'}
});
var text = ui.Label(
    'Results from analysis of Landsat images using MNDWI, EVI and NDVI',
    {fontSize: '12px'});

var startDate = ui.Textbox('YYYY-MM-DD', '2023-01-01');
var endDate =   ui.Textbox('YYYY-MM-DD', '2023-12-31');

// Define options.
var OCC = 'Water Occurrence',
    JRC = 'Transition (JRC)',
    GWL = 'Groundwater Level';

// Define the select button for options
var selection = ui.Select({
  items:[OCC,JRC,GWL],
  placeholder:'Please select an option..',
  onChange: select_f,
  value: 'Water Occurrence',
  });

mainPanel.add(selection);


//Selection function for Drop Down
function select_f() { 
  var choice = selection.getValue();
  if (choice == OCC){
      loadFreq();
  }
  else if(choice == JRC){
      loadJRC();
  }

  else if(choice == GWL){
      loadGWL();
  }
}

//=================================================   Water Occurrence  =================================================
// Define a function for displaying water occurrence
function loadFreq() {
  Map.clear();
  mainPanel.clear();
  mainPanel.add(selection);
  
  // Get the values of start and end date as client-side strings
  var startDateValue = startDate.getValue();
  var endDateValue = endDate.getValue();

  // Convert the strings to JavaScript Date objects
  var start_date = new Date(startDateValue);
  var end_date = new Date(endDateValue);
  
  // Define the valid date range as JavaScript Date objects
  var minDate = new Date('1987-01-01');
  var maxDate = new Date();

  // Check if the selected date range is out of bounds
  if (start_date < minDate || end_date > maxDate) {
    var errorMessage = ui.Label({
      value: 'Error: Date range is out of bounds. Please select a date between 1987-01-01 and current date. Relaunch app to continue.',
      style: {color: 'red', fontSize: '16px', margin: '8px 0'}
    });
    mainPanel.add(errorMessage);
    return; // Stop the function execution
  }

  // Now convert back to Earth Engine ee.Date objects for further processing
  var ee_start_date = ee.Date(startDateValue);
  var ee_end_date = ee.Date(endDateValue);

  var data_snippet, Blue, Green, Red, NIR, SWIR;

  // Compare using JavaScript Date objects
  if (start_date >= new Date('2013-03-02')) {
    data_snippet = 'LANDSAT/LC08/C02/T1_L2';
    Blue = 'SR_B2'; Green = 'SR_B3'; Red = 'SR_B4'; NIR = 'SR_B5'; SWIR = 'SR_B6';
  } else if (start_date >= new Date('2000-01-02') && start_date <= new Date('2013-03-01')) {
    data_snippet = 'LANDSAT/LE07/C02/T1_L2';
    Blue = 'SR_B1'; Green = 'SR_B2'; Red = 'SR_B3'; NIR = 'SR_B4'; SWIR = 'SR_B5';
  } else if (start_date >= new Date('1987-01-01') && start_date <= new Date('2000-01-01')) {
    data_snippet = 'LANDSAT/LT05/C02/T1_L2';
    Blue = 'SR_B1'; Green = 'SR_B2'; Red = 'SR_B3'; NIR = 'SR_B4'; SWIR = 'SR_B5';
  } else {
    var errorMessage = ui.Label({
      value: 'Error: Date range is out of bounds. Please select a valid date range.',
      style: {color: 'red', fontSize: '14px', margin: '8px 0'}
    });
    mainPanel.add(errorMessage);
    return;
  }

  // Proceed with the rest of the code for loading Landsat data
  var landsat = ee.ImageCollection(data_snippet);
  var non_mosaic = landsat.filterBounds(region).filterDate(ee_start_date, ee_end_date);

  var maskCloudsAndShadows = function(image) {
    var QA = image.select('QA_PIXEL');
    var cloudShadowBitMask = 1 << 3; 
    var cloudsBitMask = 1 << 5; 
    var cloudShadowMask = QA.bitwiseAnd(cloudShadowBitMask).eq(0);
    var cloudsMask = QA.bitwiseAnd(cloudsBitMask).eq(0);
    return image.updateMask(cloudShadowMask).updateMask(cloudsMask);
  };

  var clouds_free = non_mosaic.map(maskCloudsAndShadows);

  var water_occurrence = function(image) {
    var mndwi = image.normalizedDifference([Green, SWIR]);
    var evi = image.expression(
      '(NIR - RED)/(NIR+6*RED-7.5*BLUE+1)*2.5', {
      'BLUE': image.select(Blue),
      'RED': image.select(Red),
      'NIR': image.select(NIR),
    });
    var ndvi = image.normalizedDifference([NIR, Red]);
    var waterid = evi.lt(0.1).and(mndwi.gt(evi).or(mndwi.gt(ndvi))).rename('water');
    var water_b = waterid.clip(water_logging_prone_area);
    return water_b;
  };
  
  var result = clouds_free.map(water_occurrence);
  
  var freq = result.sum().divide(result.count());
  var msk = freq.gt(0.25); 
  var freq1 = freq.updateMask(msk);
  
  var image_w = water_poly.reduceToImage({
  properties: ['water'],
  reducer: ee.Reducer.sum()
  });
  
  var msk2 = msk.neq(image_w).and(msk.gt(0));
  
  var VisParamWar = {"bands":["water"],"min":0,"max":1, "palette": ['bee7ff','4f98f1','2652f2','091ae1']};
  // Map.centerObject(region,7.5);
  Map.setOptions('SATELLITE');
  Map.setCenter(74.00526, 29.30566, 11.5);
  Map.addLayer(freq1,VisParamWar,'Water freq');
 

//   Export.image.toDrive({
//   image: freq1,
//   description: 'Water_occurrence_frequency',
//   folder: 'Water Occurrence GEE2a',
//   region: region,
//   crs: 'EPSG:4326',
//   scale:30,
//   maxPixels:(1e12)
 
// });

  // Water-logged Area Calculation by District
  var calculateArea = function(feature) {
      var areas_img = ee.Image.pixelArea().addBands(msk2);
      var areas = areas_img.reduceRegion({
        reducer: ee.Reducer.sum()
                    .group({
                    groupField: 1,
                    groupName: 'class',
                  }),
      geometry: feature.geometry(),
      scale: 30,
      maxPixels: 1e11,
      tileScale: 2
      });
      var classAreas = ee.List(areas.get('groups'));
      var classAreaLists = classAreas.map(function(item) {
        var areaDict = ee.Dictionary(item);
        var classNumber = ee.Number(areaDict.get('class')).format();
        var area = ee.Number(areaDict.get('sum')).divide(1e6);//.round()
       return ee.List([classNumber, area]);
      });
      var result = ee.Dictionary(classAreaLists.flatten());
      // The result dictionary has area for all the classes
      // We add the district name to the dictionary and create a feature
      var district = feature.get('District');
      return ee.Feature(feature.geometry(), result.set('district', district));
  };
  var districtAreas = ignp_dist.map(calculateArea);

  var classes = ee.List.sequence(1, 1);
  // As we need to use the list of output fields in the Export function
  // we have to call .getInfo() to get the list values on the client-side
  var outputFields = ee.List(['district']).cat(classes).getInfo();

  var districtAreas_ignp = districtAreas.filter(ee.Filter.gt('1', 0));
  
  var district_list = districtAreas_ignp.reduceColumns(ee.Reducer.toList(), ['district']).get('list');
  var water_area_list = districtAreas_ignp.reduceColumns(ee.Reducer.toList(), ['1']).get('list');

  var textStyle = {
    color: 'grey',
    fontName: 'arial',
    fontSize: 8,
    rotate: 0,
    bold: false,
    italic: false
  };

  var chart = ui.Chart.array.values(ee.List(water_area_list), 0, ee.List(district_list))
              .setChartType('ColumnChart')
                .setOptions({
                    title: 'Area Water Logged',
                    intervals: {style: 'area'},
                    hAxis: {
                      title: 'District',
                      textStyle: textStyle,//{transform: 'rotate(180deg)'},
                      titleTextStyle: {italic: false, bold: true},
                    },
                    vAxis: {title: 'Area (km\u00B2)', titleTextStyle: {italic: false, bold: true}},
                    colors: ['#097b87'],
                    legend: {position: 'none'}
                });
  
  // mainPanel.widgets().set(2, chart);
  // print(chart);
  var prompt1 = ui.Label(
    'Please select the time period (from 1987 onwards)',
    {fontSize: '14px', color:'red',margin: '8px 0 -3px 8px'});

  var prompt2 = ui.Label(
    '(300+ Landsat images to be processed for each year - Longer periods can take more time to load)',
    {fontSize: '14px', color:'red'});

  var text1 = ui.Label(
    'Start date',
    { margin: '8px 0 -3px 8px',
      fontSize: '12px',
      color: 'gray'
    });

  var text2 = ui.Label(
    'End date',
    { margin: '8px 0 -3px 8px',
      fontSize: '12px',
      color: 'gray'
    });
  
  mainPanel.add(title);
  mainPanel.add(text);
  mainPanel.add(prompt1);
  mainPanel.add(prompt2);
  mainPanel.add(text1);
  mainPanel.add(startDate);
  mainPanel.add(text2);
  mainPanel.add(endDate);  
  var dropdownPanel = ui.Panel({
  layout: ui.Panel.Layout.flow('horizontal'),
  });
  mainPanel.add(dropdownPanel);
  var button2 = ui.Button('Load');
  dropdownPanel.add(button2);
  mainPanel.add(chart);
  button2.onClick(loadFreq);

  // Legend 

  // Vis parameter:  
  var vis_wff = {
    min: 25,
    max: 100,
    opacity: 1,
    palette: ['bee7ff','4f98f1','2652f2','091ae1']
  };

  // Create color bar
  function makeColorBarParams(palette) {
    return {
      bbox: [0, 0, 100, 10],
      dimensions: '120x15',
      format: 'png',
      position: 'bottom-left',
      min: 0,
      max: 100,
      palette: palette,
    };
  }

  // Thumbnail for the color bar
  var colorBar = ui.Thumbnail({
    image: ee.Image.pixelLonLat().select(0),
    params: makeColorBarParams(vis_wff.palette),
    style: {stretch: 'horizontal', margin: '0px 8px', maxHeight: '24px', position: 'bottom-left' },
  });


  // Title  
  var legendTitle = ui.Label({
    value: 'Water occurrence (%)',
    style: {fontWeight: 'bold',margin: '0px 8px', maxHeight: '24px', position:'bottom-center'}
  });
  
  
  // Labels
  var legendLabels = ui.Panel({
    widgets: [
      ui.Label(vis_wff.min, {margin: '4px 8px'}),
      ui.Label(
          ((vis_wff.max-vis_wff.min) / 2+vis_wff.min),
          {margin: '4px 8px', textAlign: 'center', stretch: 'horizontal'}),
      ui.Label(vis_wff.max, {margin: '4px 8px'})
    ],
    layout: ui.Panel.Layout.flow('horizontal')
  });


  // Add the legend to the map
  var legendPanel = ui.Panel([legendTitle, colorBar, legendLabels]);
  legendPanel.style().set({
    position: 'bottom-left'
  });
  Map.add(legendPanel);
  
}

//============================================ JRC Transition ======================================================
// Define a function for loading JRC Data
function loadJRC() {
  Map.clear();
  mainPanel.clear();
  mainPanel.add(selection);
  var dataset = ee.Image('JRC/GSW1_4/GlobalSurfaceWater');
  var jrc_water = dataset.clip(region);
    var VisParam = {"bands": ["transition"],"min":0,"max":10,"palette": 
    ['#ffffff','#0000ff','#22b14c','#d1102d','#99d9ea','#b5e61d','#e6a1aa','#ff7f27','#ffc90e','#7f7f7f','#c3c3c3']};
  
  Map.setOptions('SATELLITE');
  Map.setCenter(74.00526, 29.30566, 11.5);
  Map.addLayer(jrc_water,VisParam,'JRC Water Data');

  // Set position of panel
  var legend = ui.Panel({
  style: {
    position: 'bottom-left',
    padding: '8px 15px'
    }
  });

  // Create legend title
  var legendTitle = ui.Label({
    value: 'Transition Classes',
    style: {
    fontWeight: 'bold',
    fontSize: '18px',
    margin: '0 0 4px 0',
    padding: '0'
    }
  });

  // Add the title to the panel
  legend.add(legendTitle);
    
  // Creates and styles 1 row of the legend.
  var makeRow = function(color, name) {
      
      // Create the label that is actually the colored box.
      var colorBox = ui.Label({
        style: {
          backgroundColor:color,
          // Use padding to give the box height and width.
          padding: '8px',
          margin: '0 0 4px 0'
        }
      });
      
      // Create the label filled with the description text.
      var description = ui.Label({
        value: name,
        style: {margin: '0 0 4px 6px'}
      });
      
      // return the panel
      return ui.Panel({
        widgets: [colorBox, description],
        layout: ui.Panel.Layout.Flow('horizontal')
      });
  };


  //Palette with the colors
  var palette =['#0000ff','#22b14c','#d1102d','#99d9ea','#b5e61d','#e6a1aa','#ff7f27','#ffc90e','#7f7f7f','#c3c3c3'];

  var names = ['Permanent','New permanent','Lost permanent','Seasonal','New seasonal',
              'Lost seasonal','Seasonal to permanent','Permanent to seasonal','Ephemeral permanent','Ephemeral seasonal'];

  // Add color and names
  for (var i = 0; i < 10; i++) {
    legend.add(makeRow(palette[i], names[i]));
    }  

  // Add legend to map  
  Map.add(legend);

  // // Add TITLES AND EXPLANATORY TEXT 
  var title = ui.Label({  value: 'Water transition classes in IGNP Command',
  style: {'fontSize': '20px',color:'blue',margin: '8px 0 -3px 8px'}
  });
  var text = ui.Label(
    'These classes were generated using scenes from Landsat 5, 7, and 8 acquired between 16 March 1984 and 31 December 2021. Each pixel was individually classified into water / non-water using an expert system and the results were collated into a monthly history for the entire time period and two epochs (1984-1999, 2000-2021) for change detection. Areas where water has never been detected are masked.',
    {fontSize: '12px'});
    
  // // Create a hyperlink to an external reference.
  var link = ui.Label(
    'For more information see the associated journal article: High-resolution mapping of global surface water and its long-term changes (Nature, 2016).', {},
    'https://www.nature.com/nature/journal/v540/n7633/full/nature20584.html');
  
  mainPanel.add(title);
  mainPanel.add(text);
  mainPanel.add(link);
  
}


//=========================================   Groundwater Levels  ===========================================================
// Define a function for displaying GW Levels

function loadGWL() {
 
  Map.clear();
  Map.setCenter(72.498, 28.291, 9);
  mainPanel.clear();
  mainPanel.add(selection);
  
  var GW_bands =  ['b1','b2','b3','b4','b5','b6','b7','b8','b9','b10','b11','b12','b13','b14',
                  'b15','b16','b17','b18','b19','b20','b21','b22','b23','b24','b25','b26','b27'];
  var year = ['1994','1995','1996','1997','1998','1999','2000','2001','2002','2003','2004','2005','2006','2007','2008','2009','2010',
              '2011','2012','2013','2014','2015','2016','2017','2018','2019','2020'];
  
  //var GWL = GW_Layer.select(GW_bands).rename(year);
  var palette_gw = ['#0099ff','#ffff66','#ff3300'];
  
  
  var vis_para_gw = {
    visParams: {min: 1, max: 27, palette: palette_gw},
    defaultVisibility: true
  };
  
  // Legend 

  // Vis parameter:  
  var vis_wff = {
    min: 0,
    max: 60,
    opacity: 1,
    palette: ['#0099ff','#ffff66','#ff3300']
  };

  // Create color bar
  function makeColorBarParams(palette) {
    return {
      bbox: [0, 0, 1, 0.1],
      dimensions: '100x10',
      format: 'png',
      position: 'bottom-left',
      min: 0,
      max: 1,
      palette: palette
    };
  }

  // Thumbnail for the color bar
  var colorBar = ui.Thumbnail({
    image: ee.Image.pixelLonLat().select(0),
    params: makeColorBarParams(vis_wff.palette),
    style: {stretch: 'horizontal', margin: '0px 8px', maxHeight: '24px', position: 'bottom-left' }, 
  });


  // Title  
  var legendTitle = ui.Label({
    value: 'Groundwater Level (m-bgl)',
    style: {fontWeight: 'bold'}
  });
  
  
  // Labels
  var legendLabels = ui.Panel({
    widgets: [
      ui.Label(vis_wff.min, {margin: '4px 8px'}),
      ui.Label(
          ((vis_wff.max-vis_wff.min) / 2+vis_wff.min),
          {margin: '4px 8px', textAlign: 'center', stretch: 'horizontal'}),
      ui.Label(vis_wff.max, {margin: '4px 8px'})
    ],
    layout: ui.Panel.Layout.flow('horizontal')
  });

  // Add the legend to the map
  var legendPanel = ui.Panel([legendTitle, colorBar, legendLabels]);
  legendPanel.style().set({
    position: 'bottom-left'
  });
  

  // CREATE A SLIDER
  var selectItems = year;
  var yearSelect = ui.Slider({
    min: 1994,
    max: 2020,
    step: 1,
    style: {stretch: 'horizontal', width:'250px', position: 'top-right', whiteSpace: 'nowrap'},
      onChange: function(selected) {
        Map.clear();
        mainPanel.clear();
        mainPanel.add(selection);
        var opacitySlider = ui.Slider({
        min: 0,
        max: 1,
        value: 1,
        step: 0.01,
        });
        opacitySlider.onSlide(function(value) {
        Map.layers().forEach(function(element, index) {
        element.setOpacity(value);
        });
    });
      
        var year_n = yearSelect.getValue();
        yearSelect.setValue(year_n);
        mainPanel.add(ui.Label('Use slider to change year', {'font-size': '14px', 'fontWeight': '100'}));
        mainPanel.add(yearSelect);
        mainPanel.add(ui.Label('Use slider to change opacity', {'font-size': '14px', 'fontWeight': '100'}));
        mainPanel.add(opacitySlider);
  
        Map.add(legendPanel);
        var year = yearSelect.getValue();
        var band = ee.Number(year-1993).format('b%.0f');
        var year_str = ee.Number(year).format('%.0f').getInfo();
        var year_str_table = ee.Number(year).format('Minimum - %.0f').getInfo();
        var image_b = GW_Layer.select(band).visualize(vis_para_gw.visParams);
       
        Map.addLayer(image_b, {},year_str, vis_para_gw.defaultVisibility);
        
        var image_c = GW_Layer.select(band);
        var lines = ee.List.sequence(0, 60, 1);

        var contourlines = lines.map(function(line) {
        var mycontour = image_c
          .convolve(ee.Kernel.gaussian(5, 3))
          .subtract(ee.Image.constant(line)).zeroCrossing() 
          .multiply(ee.Image.constant(line)).toFloat();
    
        return mycontour.mask(mycontour);
        });

        contourlines = ee.ImageCollection(contourlines).mosaic();

        Map.addLayer(contourlines, {min: 0, max: 60, palette:['7d7878','0a0303']}, 'Contours'+year_str);
        
        // var chart =
        // ui.Chart.image.histogram({image: image_c, region: region, scale: 500})
        // .setSeriesNames([''])
        // .setOptions({
        //   title: 'Histogram',
        //   hAxis: {
        //     title: 'GWL (m-bgl)',
        //     titleTextStyle: {italic: false, bold: true},
        //   },
        //   vAxis:
        //       {title: 'Number of wells', titleTextStyle: {italic: false, bold: true}},
        //   colors: ['cf513e', '1d6b99', 'f0af07']
        // });
        
        // mainPanel.add(chart);
        
        var gwt2 = gwt.select('Minimum - '+year);
        
        var chart2 =
        ui.Chart.feature.histogram({features: gwt2, property: year_str_table})
        .setSeriesNames([''])
        .setOptions({
          title: 'Histogram',
          hAxis: {
            title: 'GWL (m-bgl)',
            titleTextStyle: {italic: false, bold: true},
          },
          vAxis:
              {title: 'Number of wells', titleTextStyle: {italic: false, bold: true}},
          colors: ['cf513e', '1d6b99', 'f0af07']
        });
        
        mainPanel.add(chart2);
      }
  });
  
  var opacitySlider = ui.Slider({
  min: 0,
  max: 1,
  value: 1,
  step: 0.01,
  });
  opacitySlider.onSlide(function(value) {
  Map.layers().forEach(function(element, index) {
    element.setOpacity(value);
  });
  });
  
  yearSelect.setValue(2010);
  
}
Map.setOptions('SATELLITE');
Map.setCenter(73.9917, 29.3733, 10);
ui.root.add(mainPanel);
loadFreq();
