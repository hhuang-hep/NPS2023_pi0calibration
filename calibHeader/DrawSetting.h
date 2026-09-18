void Global_setting() // global settings
{
    TGaxis::SetMaxDigits(3);
    gStyle->SetOptStat(0);
    gStyle->SetPadLeftMargin(0.13);
    gStyle->SetPadRightMargin(0.13);
    gStyle->SetLineWidth(3);
    gStyle->SetFrameLineWidth(3);
}

void format_Line(TAttLine *line, Style_t style, Color_t color, Width_t width) // Format of Line
{
    line->SetLineStyle(style);
    line->SetLineColor(color);
    line->SetLineWidth(width);
}

void format_Mark(TAttMarker *mark, Style_t style, Color_t color, Size_t msize) // Format of Marker
{
    mark->SetMarkerStyle(style);
    mark->SetMarkerColor(color);
    mark->SetMarkerSize(msize);
}

void format_Legend(TLegend *legend)  // Format of Legend
{
    legend->SetBorderSize(0);
    legend->SetFillStyle(0);
    legend->SetFillColor(0);
    legend->SetTextFont(42);
    legend->SetTextSize(0.03);
}