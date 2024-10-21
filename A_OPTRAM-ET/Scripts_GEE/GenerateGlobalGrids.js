/************************************************************
 ***        GLOBAL IRRIGATED AREA MAPPING                 ***
 ************************************************************
 *                                                          *
 *                File Name : global_grid                   *
 *               Programmer : Zhaoyuan Yao                  *
 *             Created Date : 2020-06-05                    *
 *             Last Updated : 2020-07-15                    *
 *               Contact us : yzy.sess@pku.edu.cn           *
 *                                                          *
 ************************************************************
 * File Description:                                        *
 *                                                          *
 *    The function of this script is just for dividing the  *
 *    earth into grids, which means that GEE will handle    *
 *    small areas in parallel instead of a big region that  *
 *    always results in an error: 'User Memory Limited'.    *
 *                                                          *
 *==========================================================*/
 


var boxSize   = 2;
// var lats  = ee.List.sequence(17, 19, boxSize);
// var lons  = ee.List.sequence(70, 74, boxSize);
// var lats  = ee.List.sequence(18, 20, boxSize);
// var lons  = ee.List.sequence(73, 136, boxSize);
var lats  = ee.List.sequence(-90, 90 - boxSize, boxSize);
var lons  = ee.List.sequence(-180, 180 - boxSize, boxSize);
var dataset = ee.ImageCollection('ECMWF/ERA5_LAND/DAILY_AGGR').first().select('temperature_2m');
Map.addLayer(dataset)


print(lats)
print(lons)

var grids = lats.iterate(function(lat, list){
    return ee.List(list)
             .cat(ee.List(lons.iterate(function(lon, l){
                            return ee.List(l).add(ee.Feature(ee.Geometry.Rectangle([ee.Number(lon),
                                                                                    ee.Number(lat),
                                                                                    ee.Number(lon).add(boxSize),
                                                                                    ee.Number(lat).add(boxSize)
                                                                        ])).set('id', ee.Number(ee.List(list).size()).add(ee.List(l).size()))
                                                            );
                          }
                          , ee.List([])
                          ))
              );
  }, ee.List([])); 
Map.addLayer(ee.FeatureCollection(ee.List(grids)))

// var selected = ET_dataset.reduceRegions({
//   collection: grids,
//   reducer:ee.Reducer.count(),
//   scale:500,
//   // crs, crsTransform, tileScale
// })
var selected = ee.List(grids).map(function(fea){
  var s = dataset.reduceRegion({
    reducer: ee.Reducer.count(),
    geometry:ee.Feature(fea).geometry(),
    scale: 1000, 
    bestEffort: true,
    // crs, crsTransform, maxPixels, tileScale
  })
  // return s
  return ee.Algorithms.If({
    condition: s.get('temperature_2m'),
    trueCase : fea,
    falseCase: null
  })
}, true)

print(ee.List(grids).size())
print(ee.List(grids).get(0))
print(ee.List(selected))
print(ee.List(selected).size())
Map.addLayer(ee.FeatureCollection(ee.List(selected)))

Export.table.toAsset({
  collection  : ee.FeatureCollection(ee.List(selected)),
  description : 'Global_Grid_2_degree',
  assetId     : 'Grids/Global_Grid_2_degree'
});


