import json, sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')
tpl = open('web/template.html', encoding='utf-8').read()
it = json.load(open('it.json', encoding='utf-8'))
extra = json.load(open('web/extra.json', encoding='utf-8'))
data = json.load(open('data.json', encoding='utf-8'))
tr = {}
try:
    tr = json.load(open('tr.json', encoding='utf-8'))
except FileNotFoundError:
    pass
i18n = {'it': {**tr.get('it', it), **extra['it']}, 'en': {**tr.get('en', {}), **extra['en']}, 'es': {**tr.get('es', {}), **extra['es']}}
for k in ('en', 'es'):
    missing = [x for x in i18n['it'] if x not in i18n[k]]
    print(k, 'missing', len(missing))
out = tpl.replace('/*__I18N__*/null', json.dumps(i18n, ensure_ascii=False)).replace('/*__DATA__*/null', json.dumps(data, ensure_ascii=False))
open('web/universo.html', 'w', encoding='utf-8').write(out)
print('written', len(out))
