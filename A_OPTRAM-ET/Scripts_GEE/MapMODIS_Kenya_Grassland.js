// Copyright © 2020 Zhaoyuan Yao, All right reserved.

var center = ee.Geometry.Point([34.9533, -1.6542])
var ndvi_min = 0.1;
var ndvi_max = 0.9;

var para = ee.Image('projects/ee-zhaoyuan-yao/assets/OPTRAM/MOD09A1_OPTRAM_GLOBAL_IMG_good0229_buffer200km_5')

var wet_k = para.select('top_k');
var wet_b = para.select('top_b');
var dry_k = para.select('bottom_k');
var dry_b = para.select('bottom_b');

var filterDate = ee.Filter.date('2016-06-01', '2016-07-01');
var image_ET = function(image){
  // var image = ee.ImageCollection('MODIS/061/MOD09GA').filter(filterDate).first();
  // var image = ee.ImageCollection('MODIS/061/MOD09GA').filter(filterDate).first();
  // var sr_qa = image.select('StateQA');
  var sr_qa = image.select('state_1km');
  var cloudMask = sr_qa.bitwiseAnd(1 << 0)
              .eq(sr_qa.bitwiseAnd(1 << 1)); // cloud clear
  var cloudshadowMask = sr_qa.bitwiseAnd(1 << 2).eq(0); // cloud shadow clear
  var snowiceMask = sr_qa.bitwiseAnd(1 << 12).eq(0); // snow/ice clear
  var cirrusMask  = sr_qa.bitwiseAnd(1 << 9).eq(0); // cirrus clear
  var waterMask = sr_qa.bitwiseAnd(7 << 3).eq(8); // water, bit 3-5 
  var mask = cloudMask
        .and(cloudshadowMask)
        .and(snowiceMask)
        .and(cirrusMask)
        .and(waterMask)
  var date = ee.Date(ee.Image(image).get('system:time_start'));
  // print(date, image)
  
  var red = ee.Image(image).select('sur_refl_b01').multiply(0.0001);
  var nir = ee.Image(image).select('sur_refl_b02').multiply(0.0001);
  var swir = ee.Image(image).select(['sur_refl_b07']).multiply(0.0001);
  var ndvi = nir.subtract(red).divide(nir.add(red));
  // ndvi = ndvi.multiply(ndvi.gte(ndvi_min)).add(ndvi.lt(ndvi_min).multiply(ndvi_min))
  // ndvi = ndvi.multiply(ndvi.lte(ndvi_max)).add(ndvi.gt(ndvi_max).multiply(ndvi_max))
  var fr = ndvi.subtract(ndvi_min).divide(ndvi_max - ndvi_min).pow(1);
  var str  = swir.subtract(1).multiply(swir.subtract(1)).divide(swir.multiply(2));
  
  var phi_min = fr.multiply(1.26);
  var optram = dry_b.add(dry_k.multiply(ndvi))// fr
                    .subtract(str)
                    .divide(dry_b.subtract(wet_b)
                                 .add(dry_k.subtract(wet_k).multiply(ndvi))// fr
                    )
  optram = optram.multiply(optram.gt(0))
  optram = optram.multiply(optram.lte(1)).add(optram.gt(1).multiply(1))
  var phi = optram.multiply(ee.Image(1.26).subtract(phi_min))
                  .add(phi_min);
  
  var mo_image = ee.ImageCollection("ECMWF/ERA5_LAND/DAILY_AGGR")
                   .filter(ee.Filter.date(date, date.advance(1, 'day')))
                   .first();
  var p = mo_image.select('surface_pressure').divide(1000) // Pa -> kPa
  // ACTUAL VAPOR PRESSURE [KPA]
  var tdp = mo_image.select('dewpoint_temperature_2m')
  var ea = tdp.expression(
      '0.611 * pow(10, 7.45 / (235 + Td))',{
      'Td': tdp})
  // HEAT CAPACITY OF MOIST AIT AT CONSTANT PRESSURE [J KG-1 K-1]
  var q = ee.Image(ea).multiply(0.622).divide(ee.Image(p).add(ee.Image(ea).multiply(0.622 - 1)))
  var cp = ee.Image(1).subtract(q).multiply(1003.5).add(ee.Image(q).multiply(1865))
  // SLOPE OF SATURATED VAPOR PRESWSURE VERSUS AIR TEMPERATURE [KPA K-1]
  var t_air_k = mo_image.select('temperature_2m')
  var delta = t_air_k.expression(
    '4098.0 * (0.6108 * exp(17.27 * T_C / (T_C + 237.3))) / ((T_C + 237.3) ** 2)', {
      'T_C': t_air_k.subtract(273.15)
    }).resample('bilinear')
  // LATENT HEAT OF VAPORISATION [J KG-1]
  var lambda = ee.Image(t_air_k).subtract(273.15).multiply(-2.361e-3).add(2.501).multiply(1e6)
  // PSICROMETRIC CONSTANT [KPA K-1]
  var gamma = ee.Image(cp).multiply(p).divide(ee.Image(lambda).multiply(0.622)).resample('bilinear')
  // NET RADIATION [W M-2]
  var Rn = mo_image.select('surface_net_solar_radiation_sum').add(mo_image.select('surface_net_thermal_radiation_sum')).multiply(1.1574E-5).resample('bilinear')
  
  var EF = ee.Image(phi).multiply(delta).divide(delta.add(gamma))
  
  // // v: 0.05, s: 0.4
  // var G  = fr.multiply(-1).add(1).multiply(0.4 - 0.05).add(0.05).multiply(Rn)
  // var LE = Rn.subtract(G).multiply(EF); // latent heat flux[W m-2]
  // var H  = Rn.subtract(G).multiply(EF.multiply(-1).add(1)); // sensible heat flux[W m-2]
  
  // var LE = Rn.multiply(EF).rename('LE'); // latent heat flux[W m-2]
  // var H  = Rn.multiply(EF.multiply(-1).add(1)); // sensible heat flux[W m-2]
  return Rn.multiply(EF).updateMask(mask).rename('LE')
}
var LE = ee.ImageCollection('MODIS/061/MOD09GA').filter(filterDate).map(image_ET).mean();


var visualization = {
  min: 0.0,
  max: 4,
  palette: [
    'a50026', 'd73027', 'f46d43', 'fdae61', 'fee08b', 'ffffbf',
    'd9ef8b', 'a6d96a', '66bd63', '1a9850', '006837',
  ]
};
                            
Map.addLayer(LE.divide(28.94), visualization, 'OPTRAM-ET')
Map.addLayer(ee.ImageCollection('MODIS/006/MOD16A2').select('LE').filter(filterDate).mean()
              .multiply(10000)
              .multiply(1.1574e-5)
              .divide(28.94)
              , visualization, 'MOD16')
Map.addLayer(ee.ImageCollection('CAS/IGSNRR/PML/V2_v017').filter(filterDate).mean().select('Ei')
        .add(ee.ImageCollection('CAS/IGSNRR/PML/V2_v017').filter(filterDate).mean().select('Ec'))
        .add(ee.ImageCollection('CAS/IGSNRR/PML/V2_v017').filter(filterDate).mean().select('Es'))
        , visualization, 'PML')


var export_region = ee.Geometry.BBox(32.594, -3.519, 36.395 , 0.5)
var station_name = 'Kenya';

Map.addLayer(export_region)
var ET_AOPTRAMET = LE.divide(28.94)
Export.image.toDrive({
  scale : 500,
  region: export_region,
  maxPixels: 1e13,
  crs: 'EPSG:4326',
  image: ET_AOPTRAMET,
  folder: 'global_et_example',
  description: station_name,
  fileNamePrefix: station_name + '_MOD09GA_201606',
})

var ET_MOD16A2 = ee.ImageCollection('MODIS/006/MOD16A2').select('LE').filter(filterDate).mean()
              .multiply(10000)
              .multiply(1.1574e-5)
              .divide(28.94)
Export.image.toDrive({
  scale : 500,
  region: export_region,
  maxPixels: 1e13,
  crs: 'EPSG:4326',
  image: ET_MOD16A2,
  folder: 'global_et_example',
  description: station_name,
  fileNamePrefix: station_name + '_MOD16A2_201606',
})

var ET_PMLv2 = ee.ImageCollection('CAS/IGSNRR/PML/V2_v017').filter(filterDate).mean().select('Ei')
          .add(ee.ImageCollection('CAS/IGSNRR/PML/V2_v017').filter(filterDate).mean().select('Ec'))
          .add(ee.ImageCollection('CAS/IGSNRR/PML/V2_v017').filter(filterDate).mean().select('Es'))
Export.image.toDrive({
  scale : 500,
  region: export_region,
  maxPixels: 1e13,
  crs: 'EPSG:4326',
  image: ET_PMLv2,
  folder: 'global_et_example',
  description: station_name,
  fileNamePrefix: station_name + '_PMLv2_201606',
})

var TrueColor = ee.ImageCollection('MODIS/061/MOD09GA').filter(filterDate).map(function(image){
  // var sr_qa = ee.Image(image).select('StateQA');
  var sr_qa = ee.Image(image).select('state_1km');
  var cloudMask = sr_qa.bitwiseAnd(1 << 0)
              .eq(sr_qa.bitwiseAnd(1 << 1)); // cloud clear
  var cloudshadowMask = sr_qa.bitwiseAnd(1 << 2).eq(0); // cloud shadow clear
  var snowiceMask = sr_qa.bitwiseAnd(1 << 12).eq(0); // snow/ice clear
  var cirrusMask  = sr_qa.bitwiseAnd(1 << 9).eq(0); // cirrus clear
  var waterMask = sr_qa.bitwiseAnd(7 << 3).eq(8); // water, bit 3-5 
  var mask = cloudMask
        .and(cloudshadowMask)
        // .and(snowiceMask)
        .and(cirrusMask)
        // .and(waterMask)
  var red = ee.Image(image).select('sur_refl_b01').multiply(0.0001);
  var green = ee.Image(image).select('sur_refl_b04').multiply(0.0001);
  var blue = ee.Image(image).select(['sur_refl_b03']).multiply(0.0001);
  return red.rename('red').addBands(green.rename('green')).addBands(blue.rename('blue')).updateMask(mask)
}).mean()
Export.image.toDrive({
  scale : 500,
  region: export_region,
  maxPixels: 1e13,
  crs: 'EPSG:4326',
  image: TrueColor.select(['red', 'green', 'blue']),
  folder: 'global_et_example',
  description: station_name,
  fileNamePrefix: station_name + '_TrueColor_MOD09GA_201606',
})
