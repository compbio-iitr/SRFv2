$(document).ready(function() {
    $(".redundant_button").click(function()
    {
        $(".primary_redundant").toggleClass("hide_unhide");
        $(".redundant").toggleClass("hide_unhide");	
        $(".redundant-text").toggleClass("hide_unhide");
        $(".redundant_detailed").toggleClass("hide_unhide");
        console.log('button clicked');
    });	


$(".view-winspectra").click(function()
{
    
    var image = new Image();
    var win_mer = $(this).attr('data-mer');
    var subseq = $(this).attr('data-region');
    var spectype = $(this).attr('data-spectype'); 
    var output_dir = $(this).attr('data-output');
    console.log(win_mer);
    if(win_mer >= 2)
    {
        if(spectype == "fft")
        {
            image_path = './fourier_wspec' +subseq + '_' + win_mer + '.png';
        }
        if(spectype == "dct")
        {
            image_path = './dct_wspec' +subseq + '_' + win_mer + '.png';
        }
        if(spectype == "dst")
        {
            image_path = './dst_wspec' +subseq + '_' + win_mer + '.png';
        }
        if(spectype == "strans")
        {
            image_path = './stransform_wspec' +subseq + '_' + win_mer + '.png';
        }

    }
    if(win_mer == 0) 
    {
        if(spectype == "fft")
        {
            image_path = './fourier_spectra' + subseq + '.png';
        }
        if(spectype == "dct")
        {
            image_path = './dct_spectra' + subseq + '.png';
        }
        if(spectype == "dst")
        {
            image_path = './dst_spectra' + subseq + '.png';
        }
        if(spectype == "strans")
        {
            image_path = './stransform_spectra' + subseq + '.png';
        }
   
    }
    image.src = image_path
    image.onload = function () {
        console.log('image click successful')
        var aspect = image.width/image.height;
        image.width = 600;
        image.height = Math.floor(image.width/aspect);
        $("#popup").css("width", image.width + 4);
        $("#popup").css("height", image.height + 4);
        $("#popup").css("margin-top", $(".content_container").scrollTop());
        $("#popup").css("margin-left", 50);
        $("#popup").css({ "top" : "50px", "left": "50%"}); 
        $('#popup').empty().append(image);
        $('#popup').append("<button type=\"button\" class=\"btn btn-danger btn-sm close-popup\" aria-label=\"Close\"><span aria-hidden=\"true\">X</span></button>");
        $('.close-popup').click(function()
        {
            $('#popup').fadeOut("fast", function(){});
        });
    };
    image.onerror = function ()
    {
        console.log("Error loading image");
        $('#popup').empty().append("<span class=\"loading-text\"><i class=\"glyphicon glyphicon-alert\"></i>&nbsp;Unable to load image.</span>");
        $('#popup').append("<button type=\"button\" class=\"close-popup\" aria-label=\"Close\"><span aria-hidden=\"true\">&times;</span></button>");
        $('.close-popup').click(function()
        {
            $('#popup').fadeOut("fast", function(){});
        });
    };
    $('#popup').empty().append("<span class=\"loading-text\"><i class=\"glyphicon glyphicon-refresh gly-spin\"></i>&nbsp;Loading...</span>");
    $('#popup').append("<span class=\" close-popup\"><i class=\"glyphicon glyphicon-remove\"></i></span>");
    $('.close-popup').click(function()
        {
            $('#popup').fadeOut("fast", function(){});
        });
    $('#popup').fadeIn("fast", function(){});
    });

    $(".logo_button").click(function()
    {
        var offset = $(this).offset();
        console.log(offset.top+'px');
        var image = new Image();
        var win_mer = $(this).attr('data-mer');
        var subseq = $(this).attr('data-region');
        var spectype = $(this).attr('data-spectype'); 
        var number = $(this).attr('data-number');
        var output_dir = $(this).attr('data-output');
        if(spectype == "logo") {
            image_path = './weblogo' + subseq + '_' + win_mer + '_' + number + '.png';
        }
        console.log(image_path);
        image.src = image_path; 
        image.onload = function () {
            console.log('image click successful')
            var aspect = image.width/image.height;
            image.width = Math.min(image.width, 600);
            image.height = Math.floor(image.width/aspect);
            $("#popup").css("width", image.width + 5);
            $("#popup").css("height", image.height + 5);
            $("#popup").css("position", "absolute");
            var pos_top = offset.top.toString() + 'px';
            var pos_left = offset.left.toString() + 'px';
            $("#popup").css({ "margin-top" : "-100px", "margin-left": "50px"}); 
            $("#popup").css({ "top" : pos_top, "left": pos_left}); 
            $('#popup').empty().append(image);
            $('#popup').append("<button type=\"button\" class=\"btn btn-danger btn-sm close-popup\" aria-label=\"Close\"><span aria-hidden=\"true\">X</span></button>");
            $('.close-popup').click(function()
            {
                $('#popup').fadeOut("fast", function(){});
            });
        };
        image.onerror = function ()
        {
            console.log("Error loading image");
            $('#popup').empty().append("<span class=\"loading-text\"><i class=\"glyphicon glyphicon-alert\"></i>&nbsp;Unable to load image.</span>");
            $('#popup').append("<button type=\"button\" class=\"close-popup\" aria-label=\"Close\"><span aria-hidden=\"true\">&times;</span></button>");
            $('.close-popup').click(function()
            {
                $('#popup').fadeOut("fast", function(){});
            });
        };
        $('#popup').empty().append("<span class=\"loading-text\"><i class=\"glyphicon glyphicon-refresh gly-spin\"></i>&nbsp;Loading...</span>");
        $('#popup').append("<span class=\" close-popup\"><i class=\"glyphicon glyphicon-remove\"></i></span>");
        $('.close-popup').click(function()
            {
                $('#popup').fadeOut("fast", function(){});
            });
        $('#popup').fadeIn("fast", function(){});
    });

    $('.in-div-anchor').click(function(){
        var number = $(this).attr('data-number');
        console.log('helo',number);
        $('#show-details'+number).attr("data-action", "hide");
        $('#details'+number).css("display", "block");
         var anchor = $(this).attr('href');
         var currentOffset = $('.content_container').scrollTop();
         console.log(number, 'clicked',' ',anchor,' ',$(anchor).offset().top, 'current offset',currentOffset,'value',currentOffset + $(anchor).offset().top, anchor);
         $(document).scrollTop($(anchor).offset().top);



         return false;
     });

     $('.show-details').click(function()
    {
        var region = $(this).attr('data-region');
        var action = $(this).attr('data-action');
        var div_id = "#details"+region;
        console.log("show");
        if(action == "show")
        {

            $(this).empty();
            $(this).append("<a href=\"#\"><i class=\"fa fa-angle-double-up\"></i>&nbsp;Detailed Results</a>");
            $(this).attr("data-action", "hide");
            $(div_id).css("display", "block");
        }
        else if(action == "hide")
        {
            $(this).empty();
            $(this).append("<a href=\"#\"><i class=\"fa fa-angle-double-down\"></i>&nbsp;Detailed Results</a>");
            $(this).attr("data-action", "show");
            $(div_id).css("display", "none");
        }
    });

    $(".match").each(function(i, obj)
    {
        var data = $(this).html();
        var html_string = "";
        var nucleotide_map =
        {
            'A':'adenine',
            'G':'guanine',
            'T':'thymine',
            'C':'cytosine'
        };
        var flag = 0;
        for(var i=0; i<data.length; i++)
        {
            var d = data[i];
            if(d == 'A' || d == 'C' || d == 'G' || d == 'T')
            {
                if(flag == 1)
                    html_string += "<span class=\"" + nucleotide_map[d] +"\">" + d + "</span>";
                else
                    html_string += "<span class=\"begin-pattern " + nucleotide_map[d] +"\">" + d + "</span>";
                flag = 1;
            }
            else
            {
                html_string += d;
            }
        }
        $(this).html(html_string);
    });

});