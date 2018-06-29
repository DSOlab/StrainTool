var datapoints = [];
for (i = 0; i < 500;  i++) { datapoints.push(i.toString()); }
var d_param    = 15;
var q_weights  = quadratic_weights(datapoints, d_param);
var g_weights  = gaussian_weights(datapoints, d_param);
var ctx = document.getElementById('myChart').getContext('2d');
var myChart = new Chart(ctx, {
    type: 'line',
    data: {
        labels: datapoints,
        datasets: [{
            label: 'Gaussian Weights',
            data: g_weights,
            backgroundColor: "rgba(153,255,51,0.6)"
        }, {
            label: 'Quadratic Weights',
            data: q_weights,
            backgroundColor: "rgba(255,153,0,0.6)"
        }]
    },
    options: {
        responsive: true,
        title: {
            display: true,
            text: 'Weighting For Strain'
        },
        scales: {
            xAxes: [{
                display: true,
                scaleLabel: { 
                    display: true,
                    labelString: 'Distance (km)'
                },
                ticks: {
                    maxTicksLimit: 10,
                    stepSize: 10
                }
            }],
            yAxes: [{
                display: true,
                scaleLabel: {
                    display: true,
                    labelString: 'Value'
                },
                ticks: {
                    suggestedMin: 0e0,
                    suggestedMax: 1e0,
                }
            }]
        }
    }
});
function gaussian_weights(dr, d) {
    n = dr.length;
    l = [];
    for (i = 0; i < n; i++) {
        l.push(Math.exp(-dr[i]*dr[i] / (d*d)));
    }
    return l;
}
function quadratic_weights(dr, d) {
    n = dr.length;
    l = [];
    for (i = 0; i < n; i++) {
        l.push(1e0 / (1e0 + (dr[i]*dr[i] / (d*d))));
    }
    return l;
}
function update_chart(chart, g_weights, q_weights) {
    chart.data.datasets[0].data = g_weights;
    chart.data.datasets[1].data = q_weights;
    chart.update();
}
function submit_data() {
    d_param = document.getElementById("dparam").value;
    gw = gaussian_weights(datapoints, d_param);
    qw = quadratic_weights(datapoints, d_param);
    update_chart(myChart, gw, qw);
}
/* RADAR PLOT */
var radar_datapoints = [];
var radar_data       = [];
var sta_num          = 1e0;
for (i = 0; i < 360; i++) {
    i_rad = i * Math.PI / 180e0;
    radar_data.push(sta_num*i_rad / (4e0*Math.PI));
    radar_datapoints.push(i.toString());
}
var ctx2 = document.getElementById('myChart2').getContext('2d');
var myChart2 = new Chart(ctx2, {
    type: 'line',
    data: {
        labels: radar_datapoints,
        datasets: [{
            label: 'Azimouthal Weight',
            data:radar_data,
            backgroundColor: "rgba(153,255,51,0.6)"
        }]
    },
    options: {
        responsive: true,
        title: {
            display: true,
            text: 'Weighting For Strain'
        },
        scales: {
            xAxes: [{
                display: true,
                scaleLabel: { 
                    display: true,
                    labelString: 'Theta Angle (degrees)'
                },
                ticks: {
                    maxTicksLimit: 15,
                    stepSize: 30
                }
            }],
            yAxes: [{
                display: true,
                scaleLabel: {
                    display: true,
                    labelString: 'Value'
                },
                ticks: {
                    suggestedMin: 0e0,
                    suggestedMax: 1e0,
                }
            }]
        }
    }
});
