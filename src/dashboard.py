import dash
from dash import dcc
from dash import html

import plotly.express as px
import pandas as pd
import seaborn as sns

from plotly.subplots import make_subplots
import plotly.graph_objects as go

import os

def create_dash_board():

    dir_name = os.getcwd()

    df = pd.read_csv(f"{dir_name}\Datasets\data_with_0.csv")
    df_not_null = df[df["FE, %"] > 0]
    df_raw = pd.read_csv(f"{dir_name}\\Datasets\\raw_dataset.csv")
    def sum_product():
        products = list(df["Product"].unique())
        num_products = len(products)
        colors = px.colors.qualitative.Plotly[:num_products]
        nrows = (num_products // 4) + (1 if num_products % 4 != 0 else 0)

        fig = make_subplots(rows=nrows,
                            cols=4,
                            subplot_titles=products)

        for product, c in zip(products, colors):
            act_row = products.index(product) // 4 + 1
            act_col = products.index(product) % 4 + 1
            fig.add_trace(
                go.Histogram(x=df_not_null[df_not_null["Product"] == product]["FE, %"], marker_color=c), row=act_row, col=act_col
                )

        fig.update_layout(height=600, width=1200, title_text="Распределение Фарадаевской эффективности<br>по продуктами", 
                        title_x=0.5, showlegend=False, bargap=0.2)
        return fig
    # Main metrics

    drop_data_rate = round(len(df_not_null) / len(df_raw) * 100, 2)
    null_data_rate = round(len(df_not_null) / len(df) * 100, 2)
    unic_catal = len(df.groupby(by=['DOI', 'Article name']).groups)
    sours_art = len(df_not_null["DOI"].unique())
    app = dash.Dash(__name__)

    diamonds = sns.load_dataset("diamonds")
    df = pd.read_csv(f"{dir_name}\Datasets\data_with_0.csv")

    pie = px.pie(
    data_frame=df,
    values=list(df[df["FE, %"] > 0]["Product"].value_counts()),
    names=df["Product"].unique(),
    title="Распределение ненулевых значений <br>Фарадеевской эффективности по продуктам",
    width=600,
    height=400,
    )

    pie.update_layout(
        title_x = 0.5
    )

    histogram = px.histogram(
    data_frame=df,
    x=df[df["FE, %"] > 0]["FE, %"],
    title="Общая гистограмма ненулевых значений <br>Фарадеевской эффективности",
    width=600,
    height=400,
    )

    histogram.update_layout(
        yaxis_title='FE, %',
        font=dict(size=12),
        title_x = 0.5,
        bargap=0.2
    )

    sum = sum_product()


    left_fig = html.Div(children=dcc.Graph(figure=pie))
    right_fig = html.Div(children=dcc.Graph(figure=histogram))

    top_div = html.Div(style={'backgroundColor': '#FFFFFF'}, children=[
        html.H1('Ключевые метрики проекта CO2RR', style={"font-family": "Open Sans"}),
        *[
            html.Span('Процент сохранённых после обработки строчек: ', style={"font-family": "Open Sans"}),
            html.Span(f'{drop_data_rate} %', style={'color': 'blue', "font-family": "Open Sans"}),
            html.Br(),

            html.Span('Процент нулевых строк в итоговом датасете: ', style={"font-family": "Open Sans"}),
            html.Span(f'{null_data_rate} %', style={'color': 'blue', "font-family": "Open Sans"}),
            html.Br(),

            html.Span('Количество уникальных катализаторов: ', style={"font-family": "Open Sans"}),
            html.Span(f'{unic_catal}', style={'color': 'blue', "font-family": "Open Sans"}),
            html.Br(),

            html.Span('Количество литературных источников: ', style={"font-family": "Open Sans"}),
            html.Span(f'{sours_art}', style={'color': 'blue', "font-family": "Open Sans"}),
            html.Br()          
        ]
    ])



    central_div = html.Div([left_fig, right_fig], style={"display": "flex"})
    bottom_div = html.Div(
    children=dcc.Graph(figure=sum),
    style={"display": "flex"},
    )
    app.layout = html.Div([top_div, central_div, bottom_div])

    app.run(port=8050, dev_tools_hot_reload=True)