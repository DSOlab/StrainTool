## Contributing

1. Create an issue and describe your idea
2. [Fork it](https://github.com/network#fork-destination-box)
3. Create your feature branch (`git checkout -b my-new-idea`)
4. Commit your changes (`git commit -am 'Add some feature'`)
5. Publish the branch (`git push origin my-new-idea`)
6. Create a new Pull Request
7. Profit! :white_check_mark:

### git structure

* __branch__ master/develop: stucture of main branch
  	* bin/
  		* StrainTensor.py : main script
	* data/
  		* CNRS_midas.vel : valid input file
		* station_info.dat.ref: reference output file
		* strain_info.dat: reference output file
  	* doc/
  		* 
  	* plot/
  		* gmtstrainplot.py
	* pystrain/
	* test/
	
	* .github/
  		 MarkDown templates for issues, pull requests, contributions etc.
  		


### Simple guidelines for coding

* __general__
	* Use 80 characters long line.
	* Use comments in coding.
* __gmt functions__
	* Use `-K` `-O` `-V${VRBLEVM}` at the end of each function. 
	* Add printed comments and debug messages in the code.
