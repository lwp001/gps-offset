
const X_PI:f64 = 3.14159265358979324 * 3000.0 / 180.0;
const PI:f64 = 3.1415926535897932384626;
const A:f64 = 6378245.0;
const EE:f64 = 0.00669342162296594323;

#[derive(Debug)]
pub struct Lnglat{
	//longitude
	pub lng:f64,
	//latitude
	pub lat:f64,
}


impl Lnglat{
	 pub fn new(lng:f64,lat:f64)->Self{
		Lnglat{
		  lng,
		  lat,
		}
	 }
	 //WGS84转GCj02
	 //gps to gaode(amap)
	 pub  fn wgs84_to_gcj02(&self)->Self{	 
		if out_of_china(self) {
			return Lnglat{lng:self.lng,lat:self.lat};
		}
		gcj02_wgs84_temp(self)
	 }
	 
	/**
	 * GCJ02 转换为 WGS84
	 * @param lng
	 * @param lat
	 * @returns {*[]}
	 */
	pub fn gcj02_to_wgs84(&self)->Self {
		if out_of_china(self) {
			return Lnglat{lng:self.lng,lat:self.lat};
		}
		let temp=gcj02_wgs84_temp(self);
		Lnglat{lng:self.lng * 2.0 - temp.lng,lat:self.lat * 2.0 - temp.lat}
	}
	
	
	/**
	 * 火星坐标系 (GCJ-02) 与百度坐标系 (BD-09) 的转换
	 * 即谷歌、高德 转 百度
	 */
	pub fn gcj02_to_bd09(&self)->Self {
		let z = (self.lng.powf(2.0) + self.lat.powf(2.0)).sqrt() + 0.00002 * (self.lat * X_PI).sin();
		let theta = self.lat.atan2(self.lng) + 0.000003 * (self.lng * X_PI).cos();
		Lnglat{lng: z * theta.cos() + 0.0065,
			   lat: z * theta.sin() + 0.006,
		}
	}
	
	/**
	 * 百度坐标系 (BD-09) 与 火星坐标系 (GCJ-02)的转换
	 * 即 百度 转 谷歌、高德
	 */
	pub fn bd09_to_gcj02(&self)->Self {
		let x = self.lng - 0.0065;
		let y = self.lat - 0.006;
		let z = (x * x + y * y).sqrt() - 0.00002 * (y * X_PI).sin();
		let theta = y.atan2(x) - 0.000003 * (x * X_PI).cos();
		Lnglat{lng: z * theta.cos(),
			   lat: z * theta.sin(),
		}	
	}
}

fn out_of_china(point:&Lnglat)->bool{
	return (point.lng < 72.004 || point.lng > 137.8347) || ((point.lat < 0.8293 || point.lat > 55.8271) || false);
}


fn transformlat(lng:f64, lat:f64)-> f64 {
    let mut ret:f64 = -100.0 + 2.0 * lng + 3.0 * lat + 0.2 * lat * lat + 0.1 * lng * lat + 0.2 *lng.abs().sqrt(); 
    ret += (20.0 * (6.0 * lng * PI).sin() + 20.0 * (2.0 * lng * PI).sin()) * 2.0 / 3.0;
    ret += (20.0 * (lat * PI).sin() + 40.0 * (lat / 3.0 * PI).sin()) * 2.0 / 3.0;
    ret += (160.0 *(lat / 12.0 * PI).sin() + 320.0 * (lat * PI / 30.0).sin()) * 2.0 / 3.0;
    return ret;
}
 
fn transformlng(lng:f64, lat:f64)->f64 {
    let mut ret:f64 = 300.0 + lng + 2.0 * lat + 0.1 * lng * lng + 0.1 * lng * lat + 0.1 * lng.abs().sqrt();
    ret += (20.0 * (6.0 * lng * PI).sin() + 20.0 * (2.0 * lng * PI).sin()) * 2.0 / 3.0;
    ret += (20.0 * (lng * PI).sin() + 40.0 * (lng / 3.0 * PI).sin()) * 2.0 / 3.0;
    ret += (150.0 *(lng / 12.0 * PI).sin() + 300.0 * (lng / 30.0 * PI).sin()) * 2.0 / 3.0;
    return ret;	
}


fn gcj02_wgs84_temp(point:&Lnglat)->Lnglat{
		let mut dlat:f64 = transformlat(point.lng - 105.0, point.lat - 35.0);
		let mut dlng:f64 = transformlng(point.lng - 105.0, point.lat - 35.0);
		let radlat:f64 = point.lat / 180.0 * PI;
		let mut magic:f64 = radlat.sin();
		magic = 1.0 - EE * magic * magic;
		let sqrtmagic = magic.sqrt();
		dlat = (dlat * 180.0) / ((A * (1.0 - EE)) / (magic * sqrtmagic) * PI);
		dlng = (dlng * 180.0) / (A / sqrtmagic * radlat.cos() * PI);
		Lnglat{lat: point.lat + dlat,
			   lng: point.lng + dlng,
			   }
}

