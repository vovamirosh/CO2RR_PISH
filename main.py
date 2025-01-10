from src import preprocessing, descriptors_add, sql_converting, dashboard

# Скачивание датасета из гугл таблицы, его обработкаб сохранение в обработанный файл
link_to_project = 'https://docs.google.com/spreadsheets/d/e/2PACX-1vQxIFJUE5JeauvxO11xI_LXezw_LNB_Rv4CNYoDJ0EIFKkLJiZ4ERt4CU4V5R2mEQvCp-n_A2OqKwx_/pub?gid=1053270247&single=true&output=csv'
preprocessing.preprocess(link_to_project)

# Добавление дескрипторов для неорганических соединений в датасете, сохранение нового датасета, в том числе
# расширенного датасета с 0 значениями эффективности катлизатор
descriptors_add.descr_adding()

# Создание базы данных формата sqlite
sql_converting.create_sq()

# Создание дашборда проекта. По умолчанию запускается по ссылке http://127.0.0.1:8050
dashboard.create_dash_board()